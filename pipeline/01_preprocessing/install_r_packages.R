# Install all R packages needed for the full pipeline.
# Run once before any other R scripts:
#   Rscript pipeline/01_preprocessing/install_r_packages.R

cat("Installing CRAN packages...\n")
cran_packages <- c(
  "Seurat",       # core single-cell analysis
  "ggplot2",      # plotting
  "dplyr",        # data manipulation
  "patchwork",    # combining figures
  "Matrix",       # sparse matrices
  "scales",       # axis formatting
  "viridis",      # colour scales
  "RColorBrewer", # colour palettes
  "ggrepel",      # non-overlapping labels (Phase 4/5)
  "tidyr"         # unnest() used in Phase 3
)
install.packages(cran_packages, repos = "https://cloud.r-project.org", quiet = TRUE)

cat("Installing BiocManager + Bioconductor packages...\n")
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager", repos = "https://cloud.r-project.org")

bioc_packages <- c(
  "BiocGenerics",
  "GenomicRanges",   # needed by InferCNV
  "infercnv",        # CNV calling (Phase 2)
  "biomaRt"          # Ensembl paralog query (Phase 5)
)
BiocManager::install(bioc_packages, ask = FALSE, update = FALSE)

cat("Installing GitHub packages...\n")
if (!requireNamespace("remotes", quietly = TRUE))
  install.packages("remotes", repos = "https://cloud.r-project.org")

remotes::install_github("navinlabcode/copykat",    upgrade = "never", quiet = TRUE)
remotes::install_github("navinlabcode/CopyKit",    upgrade = "never", quiet = TRUE)

cat("\nAll packages installed. Verifying...\n")
required <- c("Seurat", "ggplot2", "dplyr", "patchwork", "infercnv")
for (pkg in required) {
  ok <- requireNamespace(pkg, quietly = TRUE)
  cat(sprintf("  %-15s %s\n", pkg, if (ok) "OK" else "MISSING"))
}
