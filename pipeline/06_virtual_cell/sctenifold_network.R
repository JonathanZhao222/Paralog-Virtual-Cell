# Phase 6a — scTenifold Gene Regulatory Network
#
# Build a gene regulatory network (GRN) from the scRNA-seq data using scTenifold.
# scTenifold constructs a per-cell tensor of coexpression, then aggregates to a
# single network via tensor decomposition.
#
# We build TWO networks:
#   1. Arm-loss cells  (consensus aneuploid tumour cells with ≥1 lost arm)
#   2. Normal-reference cells
#
# The networks are saved for use in Phase 6b (virtual KO).
#
# Input:
#   data/processed/seurat_qc.rds
#   results/consensus/cell_cnv_status.csv
#   results/paralog_pairs/cross_arm_paralogs.csv
#
# Output:
#   results/ko_scores/grn_arm_loss.rds     — GRN for arm-loss cells
#   results/ko_scores/grn_normal.rds       — GRN for normal cells
#   results/ko_scores/gene_list.rds        — shared gene universe

library(Seurat)
library(scTenifoldNet)   # installs as scTenifoldNet from CRAN
library(dplyr)

# ── Paths ──────────────────────────────────────────────────────────────────────
ROOT        <- getwd()
SEURAT_RDS  <- file.path(ROOT, "data", "processed", "seurat_qc.rds")
CNV_STATUS  <- file.path(ROOT, "results", "consensus", "cell_cnv_status.csv")
PARALOGS    <- file.path(ROOT, "results", "paralog_pairs", "cross_arm_paralogs.csv")
OUT_DIR     <- file.path(ROOT, "results", "ko_scores")
dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)

# Gene universe: limit to the paralog pairs' genes + their paralogs to keep
# network tractable. scTenifold is memory-intensive for > ~2,000 genes.
MAX_GENES    <- 2000
N_CELLS_NET  <- 500    # subsample to keep tensor decomposition feasible
set.seed(42)

cat("=== Phase 6a: scTenifold Gene Regulatory Networks ===\n\n")

# ── 1. Load data ───────────────────────────────────────────────────────────────
cat("[1/4] Loading data...\n")
seurat      <- readRDS(SEURAT_RDS)
cell_status <- read.csv(CNV_STATUS)
paralogs    <- read.csv(PARALOGS)

# Cell sets
loss_cells   <- cell_status$cell[cell_status$consensus_aneuploid]
normal_cells <- colnames(seurat)[seurat$tumour_status == "normal_reference"]
loss_cells   <- intersect(loss_cells,   colnames(seurat))
normal_cells <- intersect(normal_cells, colnames(seurat))

cat(sprintf("      Arm-loss cells: %d | Normal cells: %d\n",
            length(loss_cells), length(normal_cells)))

# Subsample if needed
if (length(loss_cells)   > N_CELLS_NET) loss_cells   <- sample(loss_cells,   N_CELLS_NET)
if (length(normal_cells) > N_CELLS_NET) normal_cells <- sample(normal_cells, N_CELLS_NET)

# ── 2. Define gene universe ────────────────────────────────────────────────────
cat("[2/4] Defining gene universe...\n")

# The gene universe must contain BOTH the lost-arm gene AND its paralog for the
# KO simulation to work. Naive union fills the budget with lost-arm genes only.
# Fix: load Phase 5b upregulation results and select top significant pairs,
# ensuring both members of each pair enter the universe together.
UPREG_PATH <- file.path(ROOT, "results", "paralog_pairs", "paralog_upregulation.csv")
upreg <- read.csv(UPREG_PATH)
upreg <- upreg[order(upreg$padj, -upreg$lfc_paralog), ]
sig_pairs <- upreg[!is.na(upreg$padj) & upreg$padj < 0.05 & upreg$lfc_paralog > 0.3, ]

# Keep only pairs where both genes are measured in the Seurat object
measurable <- sig_pairs[
  sig_pairs$gene        %in% rownames(seurat) &
  sig_pairs$paralog_gene %in% rownames(seurat), ]

# Select top pairs to fill ~half the gene budget (leaving room for HVGs)
n_pairs    <- min(floor(MAX_GENES / 2), nrow(measurable))
top_pairs  <- head(measurable, n_pairs)
pair_genes <- unique(c(top_pairs$gene, top_pairs$paralog_gene))
cat(sprintf("      Top significant pairs: %d  → %d unique genes\n",
            n_pairs, length(pair_genes)))

# Fill remaining budget with highly variable genes
seurat <- FindVariableFeatures(seurat, selection.method = "vst",
                               nfeatures = MAX_GENES, verbose = FALSE)
hvg <- VariableFeatures(seurat)

gene_universe <- unique(c(pair_genes, hvg))
gene_universe <- intersect(gene_universe, rownames(seurat))
gene_universe <- head(gene_universe, MAX_GENES)
cat(sprintf("      Gene universe: %d genes (%d paired, %d HVG fill)\n",
            length(gene_universe),
            sum(gene_universe %in% pair_genes),
            sum(!gene_universe %in% pair_genes)))

saveRDS(gene_universe, file.path(OUT_DIR, "gene_list.rds"))

# ── 3. Extract count matrices ─────────────────────────────────────────────────
cat("[3/4] Extracting count matrices...\n")

raw_counts <- GetAssayData(seurat, assay = "RNA", layer = "counts")

mat_loss   <- as.matrix(raw_counts[gene_universe, loss_cells])
mat_normal <- as.matrix(raw_counts[gene_universe, normal_cells])

cat(sprintf("      Arm-loss matrix:  %d genes × %d cells\n",
            nrow(mat_loss), ncol(mat_loss)))
cat(sprintf("      Normal matrix:    %d genes × %d cells\n",
            nrow(mat_normal), ncol(mat_normal)))

# ── 4. Build GRNs with scTenifoldNet ──────────────────────────────────────────
cat("[4/4] Building GRNs (this can take 30-90 min per network)...\n")
cat("      Building arm-loss network...\n")

grn_loss <- scTenifoldNet(
  X             = mat_loss,
  Y             = mat_loss,
  nc_nNet       = 10,
  nc_nCells     = 100,
  nc_scale      = TRUE,
  td_K          = 3,
  qc_minLibSize = 1
)
saveRDS(grn_loss, file.path(OUT_DIR, "grn_arm_loss.rds"))
cat("      Saved → results/ko_scores/grn_arm_loss.rds\n")

cat("      Building normal-reference network...\n")
grn_normal <- scTenifoldNet(
  X             = mat_normal,
  Y             = mat_normal,
  nc_nNet       = 10,
  nc_nCells     = 100,
  nc_scale      = TRUE,
  td_K          = 3,
  qc_minLibSize = 1
)
saveRDS(grn_normal, file.path(OUT_DIR, "grn_normal.rds"))
cat("      Saved → results/ko_scores/grn_normal.rds\n")

cat("\n[done] Phase 6a complete.\n")
cat("Next: Rscript pipeline/06_virtual_cell/ko_simulation.R\n")
