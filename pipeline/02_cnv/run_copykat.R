# Phase 2b — CopyKAT
#
# CopyKAT detects aneuploidy by comparing local gene expression to a smooth
# baseline, using diploid (normal) cells as a reference.
#
# Input:   data/processed/seurat_qc.rds
# Output:  results/cnv_calls/copykat/
#          figures/02_cnv/copykat_heatmap.pdf

library(Seurat)
library(copykat)
library(ggplot2)

# ── Paths ─────────────────────────────────────────────────────────────────────
ROOT       <- getwd()
SEURAT_RDS <- file.path(ROOT, "data", "processed", "seurat_qc.rds")
OUT_DIR    <- file.path(ROOT, "results", "cnv_calls", "copykat")
FIG_DIR    <- file.path(ROOT, "figures", "02_cnv")
dir.create(OUT_DIR,  showWarnings = FALSE, recursive = TRUE)
dir.create(FIG_DIR,  showWarnings = FALSE, recursive = TRUE)

# ── 1. Load Seurat object ─────────────────────────────────────────────────────
cat("[1/3] Loading Seurat object...\n")
seurat <- readRDS(SEURAT_RDS)
cat(sprintf("      %d cells loaded\n", ncol(seurat)))

# ── 2. Run CopyKAT (skip if results already exist) ───────────────────────────
pred_path <- file.path(OUT_DIR, "breast_cancer_copykat_prediction.txt")

if (file.exists(pred_path)) {
  cat("[2/3] Existing CopyKAT results found — skipping computation.\n")
  cat("      Reading predictions from saved file...\n")
  pred <- read.table(pred_path, header = TRUE, sep = "\t")
  colnames(pred) <- c("cell.names", "copykat.pred")
} else {
  cat("[2/3] Running CopyKAT (2-5 hrs)...\n")
  raw_counts   <- GetAssayData(seurat, assay = "RNA", layer = "counts")
  raw_mat      <- as.matrix(raw_counts)
  normal_cells <- colnames(seurat)[seurat$tumour_status == "normal_reference"]
  cat(sprintf("      %d total cells | %d normal reference cells\n",
              ncol(seurat), length(normal_cells)))

  setwd(OUT_DIR)
  copykat_result <- copykat(
    rawmat          = raw_mat,
    id.type         = "S",
    cell.line       = "no",
    ngene.chr       = 5,
    win.size        = 25,
    KS.cut          = 0.1,
    sam.name        = "breast_cancer",
    distance        = "euclidean",
    norm.cell.names = normal_cells,
    output.seg      = FALSE,
    plot.genes      = TRUE,
    genome          = "hg20",
    n.cores         = 4
  )
  setwd(ROOT)

  pred <- if (!is.null(copykat_result$prediction)) {
    copykat_result$prediction
  } else {
    read.table(pred_path, header = TRUE, sep = "\t",
               col.names = c("cell.names", "copykat.pred"))
  }
}

# ── 3. Generate figures and save ──────────────────────────────────────────────
cat("[3/3] Processing results...\n")

# Summarise aneuploid calls per cell (used in consensus step)
aneuploid_calls <- data.frame(
  cell         = pred$cell.names,
  copykat_call = pred$copykat.pred,
  row.names    = NULL
)
write.csv(aneuploid_calls,
          file.path(OUT_DIR, "copykat_aneuploid_calls.csv"),
          row.names = FALSE)

# Figure: annotate UMAP with CopyKAT ploidy calls
seurat$copykat_ploidy <- aneuploid_calls$copykat_call[
  match(colnames(seurat), aneuploid_calls$cell)
]
seurat$copykat_ploidy[is.na(seurat$copykat_ploidy)] <- "unassigned"

p <- DimPlot(seurat, reduction = "umap", group.by = "copykat_ploidy",
             cols = c("aneuploid"  = "#E74C3C",
                      "diploid"    = "#3498DB",
                      "unassigned" = "grey80")) +
     ggtitle("CopyKAT: Aneuploid vs Diploid Calls") +
     theme_classic()

ggsave(file.path(FIG_DIR, "copykat_umap.pdf"), p, width = 7, height = 6)

n_aneuploid <- sum(aneuploid_calls$copykat_call == "aneuploid", na.rm = TRUE)
n_diploid   <- sum(aneuploid_calls$copykat_call == "diploid",   na.rm = TRUE)
cat(sprintf("      Aneuploid: %d | Diploid: %d\n", n_aneuploid, n_diploid))
cat(sprintf("[done] Results saved → results/cnv_calls/copykat/\n"))
cat("\nNext: Rscript pipeline/02_cnv/run_copykit.R\n")
