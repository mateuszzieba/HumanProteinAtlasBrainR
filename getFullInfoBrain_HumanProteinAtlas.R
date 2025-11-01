getFullInfoBrain_HumanProteinAtlas <- function(
    gene_symbol,
    section = "brain",
    retries = 5,
    wait_sec = 5,
    save_file = TRUE,
    verbose = TRUE
) {
  # --- Check required packages ---
  required_pkgs <- c(
    "httr", "jsonlite", "readr", "rvest",
    "dplyr", "stringr", "glue", "purrr", "tibble"
  )
  
  missing_pkgs <- required_pkgs[!sapply(required_pkgs, requireNamespace, quietly = TRUE)]
  if (length(missing_pkgs) > 0) {
    message("❌ Missing required packages: ", paste(missing_pkgs, collapse = ", "))
    message("Install them with:\ninstall.packages(c('", 
            paste(missing_pkgs, collapse = "', '"), "'))")
    return(NULL)
  }
  
  suppressPackageStartupMessages({
    lapply(required_pkgs, library, character.only = TRUE)
  })
  
  # --- Safe try helper ---
  safe_try <- function(expr) {
    tryCatch(expr, error = function(e) { 
      if (verbose) message("⚠️ Error: ", e$message)
      return(NULL)
    })
  }
  
  # --- Retry helper for HTTP requests ---
  retry_request <- function(url, retries, wait_sec) {
    for (i in seq_len(retries)) {
      res <- try(httr::GET(url), silent = TRUE)
      if (!inherits(res, "try-error") && httr::status_code(res) == 200) {
        return(res)
      }
      message(glue::glue("⏳ Attempt {i}/{retries} failed, waiting {wait_sec}s before retry..."))
      Sys.sleep(wait_sec)
    }
    stop(glue::glue("❌ Failed to fetch after {retries} attempts."))
  }
  
  # --- Step 1: Fetch HTML page from HPA ---
  html_text <- safe_try({
    search_url <- sprintf("https://www.proteinatlas.org/search/%s?format=json",
                          URLencode(gene_symbol, reserved = TRUE))
    res <- retry_request(search_url, retries, wait_sec)
    js <- jsonlite::fromJSON(rawToChar(res$content))
    if (length(js) == 0) stop("Gene not found: ", gene_symbol)
    df <- as.data.frame(js)
    ensg_col <- intersect(names(df), c("Ensembl","ensembl"))[1]
    ensg <- df[[ensg_col]][1]
    url <- sprintf("https://www.proteinatlas.org/%s-%s/%s", ensg, gene_symbol, section)
    message("➡️ Fetching: ", url)
    html_resp <- retry_request(url, retries, wait_sec)
    httr::content(html_resp, "text", encoding = "UTF-8")
  })
  if (is.null(html_text)) return(NULL)
  
  if (save_file) {
    file_name <- sprintf("HPA_%s_%s.html", gene_symbol, section)
    readr::write_file(html_text, file_name)
    if (verbose) message("💾 Saved: ", file_name)
  }
  
  page <- safe_try(rvest::read_html(html_text))
  if (is.null(page)) return(NULL)
  
  # --- Subfunctions ---
  extract_hpa_human_brain_dataset <- function(page) safe_try({
    span_node <- page %>% rvest::html_elements(xpath = "//span[contains(., 'HPA Human brain dataset')]")
    if (length(span_node) == 0) return(NULL)
    parent_node <- span_node %>% rvest::html_element(xpath = "./ancestor::table[1]")
    nodes <- parent_node %>% rvest::html_nodes("a.brainrna")
    if (length(nodes) == 0) return(NULL)
    tibble::tibble(
      region = nodes %>% rvest::html_text(trim = TRUE),
      brainregion_id = nodes %>% rvest::html_attr("brainregion"),
      nTPM = as.numeric(nodes %>% rvest::html_attr("nx")),
      color = nodes %>% rvest::html_attr("color")
    ) %>% dplyr::filter(!is.na(nTPM)) %>% dplyr::arrange(dplyr::desc(nTPM))
  })
  
  extract_json_dataset <- function(page, chart_id, value_name) safe_try({
    script_text <- page %>%
      rvest::html_elements("script") %>%
      rvest::html_text2() %>%
      stringr::str_subset(chart_id) %>%
      dplyr::first()
    if (is.na(script_text)) return(NULL)
    json_match <- stringr::str_match(script_text, "barChart\\((\\[\\{.*?\\}\\])")[, 2]
    if (is.na(json_match)) return(NULL)
    df <- jsonlite::fromJSON(json_match) %>% tibble::as_tibble()
    if (value_name == "enrichment_change") {
      df %>% dplyr::transmute(cell_type = label,
                              enrichment_change = as.numeric(value),
                              color) %>% dplyr::arrange(dplyr::desc(enrichment_change))
    } else {
      df %>% dplyr::transmute(
        region = label,
        brainregion_id = URLdecode(stringr::str_extract(url, "(?<=tissue\\/)[^#]+")),
        !!value_name := as.numeric(value),
        color
      ) %>% dplyr::arrange(dplyr::desc(!!rlang::sym(value_name)))
    }
  })
  
  clean_hpa_brain_table <- function(tbl_raw) safe_try({
    names(tbl_raw) <- make.names(names(tbl_raw), unique = TRUE)
    tbl <- tbl_raw %>% dplyr::select(where(~ !all(is.na(.x) | .x == "")))
    cluster_row <- tbl %>%
      dplyr::filter(dplyr::if_any(dplyr::everything(), ~ stringr::str_detect(.x, "is part of"))) %>%
      dplyr::slice_head(n = 1)
    if (nrow(cluster_row) == 0) return(NULL)
    txt <- paste(cluster_row, collapse = " ")
    cluster_info <- tibble::tibble(
      gene = stringr::str_extract(txt, "^[A-Z0-9]+"),
      cluster_number = stringr::str_extract(txt, "(?<=cluster )\\d+"),
      cluster_description = stringr::str_match(txt, "cluster \\d+ (.*?) with")[,2],
      confidence = stringr::str_extract(txt, "\\d+\\.\\d+") %>% as.numeric(),
      n_genes_in_cluster = stringr::str_extract(txt, "\\d+(?= genes in cluster)") %>% as.numeric()
    )
    correlated_genes <- tbl_raw %>%
      .[-c(1:3), c(1:4)] %>%
      as.data.frame() %>%
      setNames(c("hgnc_symbol", "Description", "Correlation", "Cluster")) %>%
      dplyr::mutate(
        Correlation = suppressWarnings(as.numeric(Correlation)),
        Cluster = suppressWarnings(as.numeric(Cluster))
      ) %>%
      dplyr::filter(!is.na(hgnc_symbol) & hgnc_symbol != "")
    list(cluster_info = cluster_info, correlated_genes = correlated_genes)
  })
  
  clean_tbl3_fn <- function(tbl_raw) safe_try({
    tbl_raw %>%
      dplyr::rename(col1 = 1, col2 = 2) %>%
      dplyr::mutate(
        across(everything(), ~ stringr::str_replace_all(., "[\\t\\n]+", " ")),
        col1 = stringr::str_replace(col1, "i\\s.*", ""),
        col1 = stringr::str_trim(col1),
        col2 = stringr::str_trim(col2)
      ) %>%
      dplyr::rename(attribute = col1, value = col2)
  })
  
  clean_tbl4_fn <- function(tbl_raw) safe_try({
    tbl_raw %>%
      tibble::as_tibble(.name_repair = "unique") %>%
      dplyr::select(1:4) %>%
      dplyr::rename_with(~ c("col1", "col2", "col3", "col4")) %>%
      dplyr::mutate(across(everything(), ~ stringr::str_replace_all(., "[\\t\\n]+", " "))) %>%
      dplyr::mutate(col1 = stringr::str_replace(col1, "i\\s.*", "")) %>%
      dplyr::mutate(across(everything(), stringr::str_trim)) %>%
      dplyr::filter(!col1 %in% c("")) %>%
      dplyr::rename(attribute = col1, human = col2, pig = col3, mouse = col4)
  })
  
  # --- Step 2: Extract datasets ---
  tables <- safe_try(page %>% rvest::html_nodes("table"))
  
  GENERAL_INFORMATION <- NULL
  HUMAN_PROTEIN_ATLAS_INFORMATION <- NULL
  
  if (!is.null(tables) && length(tables) >= 4) {
    tbl3_raw <- safe_try(rvest::html_table(tables[[3]], header = TRUE, fill = TRUE))
    tbl4_raw <- safe_try(rvest::html_table(tables[[4]], header = TRUE, fill = TRUE))
    GENERAL_INFORMATION <- clean_tbl3_fn(tbl3_raw)
    HUMAN_PROTEIN_ATLAS_INFORMATION <- clean_tbl4_fn(tbl4_raw)
  }
  
  correlated <- NULL
  if (!is.null(tables) && length(tables) >= 12) {
    tbl_raw <- safe_try(rvest::html_table(tables[[12]], header = TRUE, fill = TRUE))
    correlated <- clean_hpa_brain_table(tbl_raw)
  }
  
  # --- Step 3: Build final structured list ---
  results <- list(
    gene_symbol = gene_symbol,
    GENERAL_INFORMATION = GENERAL_INFORMATION,
    HUMAN_PROTEIN_ATLAS_INFORMATION = HUMAN_PROTEIN_ATLAS_INFORMATION,
    human_brain = extract_hpa_human_brain_dataset(page),
    pig_brain = extract_json_dataset(page, "RNAChart66", "nTPM"),
    mouse_brain = extract_json_dataset(page, "RNAChart46", "nTPM"),
    stereo_seq = extract_json_dataset(page, "RNAChart105", "enrichment_change"),
    gtex = extract_json_dataset(page, "RNAChart68", "nTPM"),
    fantom5 = extract_json_dataset(page, "RNAChart69", "scaled_TPM"),
    allen_ish = extract_json_dataset(page, "RNAChart67", "expression_energy"),
    correlated = correlated
  )
  
  if (verbose) message(glue::glue("🎯 Completed processing for {gene_symbol}."))
  return(results)
}
