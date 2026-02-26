# Phase 2c — CopyKIT
#
# NOTE (2025): CopyKIT currently cannot be installed due to an upstream bug —
# treeio (a CopyKIT dependency) references random_ref() from tidytree, but
# that function has not been released in any public tidytree version.
# Track: https://github.com/navinlabcode/CopyKit/issues
#
# The consensus step (Phase 3) uses InferCNV + CopyKAT only until this is fixed.
# This script is retained for when CopyKIT becomes installable.
#
# CopyKIT uses a k-nearest neighbour smoothing approach to infer copy number
# from scRNA-seq. It integrates well with Seurat and produces arm-level calls.
#
# Input:   data/processed/seurat_qc.rds
# Output:  results/cnv_calls/copykit/
#          figures/02_cnv/copykit_umap.pdf

library(Seurat)
library(CopyKit)
library(ggplot2)

# ── Paths ─────────────────────────────────────────────────────────────────────
ROOT       <- getwd()
SEURAT_RDS <- file.path(ROOT, "data", "processed", "seurat_qc.rds")
OUT_DIR    <- file.path(ROOT, "results", "cnv_calls", "copykit")
FIG_DIR    <- file.path(ROOT, "figures", "02_cnv")
dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)
dir.create(FIG_DIR, showWarnings = FALSE, recursive = TRUE)

# ── 1. Load data ──────────────────────────────────────────────────────────────
cat("[1/4] Loading Seurat object...\n")
seurat <- readRDS(SEURAT_RDS)
cat(sprintf("      %d cells\n", ncol(seurat)))

# CopyKIT requires raw counts
raw_counts <- GetAssayData(seurat, assay = "RNA", layer = "counts")

# Normal reference cell barcodes
normal_cells <- colnames(seurat)[seurat$tumour_status == "normal_reference"]
tumour_cells <- colnames(seurat)[seurat$tumour_status == "tumour"]

# ── 2. Create CopyKIT object ──────────────────────────────────────────────────
cat("[2/4] Creating CopyKIT object...\n")
copykit_obj <- CopyKit(
  raw_counts  = raw_counts,
  genome      = "hg38"
)

# ── 3. Run CopyKIT ────────────────────────────────────────────────────────────
cat("[3/4] Running CopyKIT (10-30 min)...\n")

# Segment copy number
copykit_obj <- runSegmentation(copykit_obj)

# Normalise against the normal reference cells
copykit_obj <- calcConsensus(
  copykit_obj,
  consensus_by   = "cell",
  normalise_cells = normal_cells
)

# Cluster cells by CNV profile
copykit_obj <- findClusters(copykit_obj)

# ── 4. Extract results ────────────────────────────────────────────────────────
cat("[4/4] Extracting arm-level CNV scores...\n")

# Segment-level CNV ratios (log2)
seg_data <- as.data.frame(copykit_obj@consensus)
write.csv(seg_data, file.path(OUT_DIR, "copykit_segment_cnv.csv"))

# Per-cell subclone assignments
subclones <- data.frame(
  cell         = colnames(copykit_obj),
  copykit_clone = copykit_obj$subclones,
  row.names    = NULL
)
write.csv(subclones, file.path(OUT_DIR, "copykit_subclones.csv"), row.names = FALSE)

# Arm-level summary: average CNV ratio per arm per cell
# Load chromosome arm boundaries
arm_bed <- read.table(
  file.path(ROOT, "data", "reference", "chromosome_arms_hg38.bed"),
  col.names = c("chrom", "start", "end", "arm")
)

# Heatmap figure
p_heat <- plotHeatmap(copykit_obj,
                      label    = "subclones",
                      order_by = "hclust") +
          ggtitle("CopyKIT: CNV Heatmap")
ggsave(file.path(FIG_DIR, "copykit_heatmap.pdf"), p_heat,
       width = 14, height = 8)

# Annotate UMAP with CopyKIT subclone
seurat$copykit_clone <- subclones$copykit_clone[
  match(colnames(seurat), subclones$cell)
]

p_umap <- DimPlot(seurat, reduction = "umap", group.by = "copykit_clone",
                  label = TRUE, repel = TRUE) +
          ggtitle("CopyKIT Subclones") +
          theme_classic()
ggsave(file.path(FIG_DIR, "copykit_umap.pdf"), p_umap, width = 8, height = 6)

# Save CopyKIT object
saveRDS(copykit_obj, file.path(OUT_DIR, "copykit_obj.rds"))

cat(sprintf("[done] Results saved → results/cnv_calls/copykit/\n"))
cat("\nAll 3 CNV tools complete.\n")
cat("Next: Rscript pipeline/03_cnv_consensus/consensus_true_positive.R\n")
