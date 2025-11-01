# 🧠 HumanProteinAtlasBrainR

**HumanProteinAtlasBrainR** is a lightweight R utility that automates extraction of *brain-specific datasets* for any gene directly from the [Human Protein Atlas (HPA)](https://www.proteinatlas.org/).

It retrieves:
- General gene information and species-specific expression tables (human, pig, mouse)
- RNA expression datasets (HPA, GTEx, FANTOM5, Allen ISH)
- Cluster-correlated genes
- All safely wrapped in one call, with retry logic and automatic package checking

---

## 🔧 Quick start (direct from GitHub)

You can load the function directly into R without downloading the repo manually:

```r
# Install 'remotes' package if not yet installed
install.packages("remotes")

# Load the function directly from GitHub (main branch)
source("https://raw.githubusercontent.com/<YOUR_USERNAME>/HumanProteinAtlasBrainR/main/getFullInfoBrain_HumanProteinAtlas.R")

# Example usage:
res <- getFullInfoBrain_HumanProteinAtlas("TSC22D3")
