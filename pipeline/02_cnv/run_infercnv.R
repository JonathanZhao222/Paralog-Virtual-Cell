# Phase 2a — InferCNV
#
# InferCNV infers CNV by comparing tumour cell expression to a normal
# reference (immune/stromal cells) across genomic windows.
#
# Input:   data/processed/seurat_qc.rds
#          data/reference/gene_order_hg38.txt
# Output:  results/cnv_calls/infercnv/    (InferCNV outputs + HMM states)
#          figures/02_cnv/infercnv_heatmap.pdf

library(Seurat)
library(infercnv)
library(ggplot2)

# ── Paths ─────────────────────────────────────────────────────────────────────
ROOT       <- getwd()
SEURAT_RDS <- file.path(ROOT, "data", "processed", "seurat_qc.rds")
GENE_ORDER <- file.path(ROOT, "data", "reference", "gene_order_hg38.txt")
OUT_DIR    <- file.path(ROOT, "results", "cnv_calls", "infercnv")
FIG_DIR    <- file.path(ROOT, "figures", "02_cnv")
dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)
dir.create(FIG_DIR, showWarnings = FALSE, recursive = TRUE)

# ── 1. Load Seurat object ─────────────────────────────────────────────────────
cat("[1/4] Loading Seurat object...\n")
seurat <- readRDS(SEURAT_RDS)
cat(sprintf("      %d cells loaded\n", ncol(seurat)))

# ── Subsample to reduce memory usage ──────────────────────────────────────────
# InferCNV is memory-intensive. Use all normal reference cells (needed for
# accurate baseline) and subsample tumour cells to a manageable number.
set.seed(42)
MAX_TUMOUR  <- 3000   # subsample tumour cells
MAX_NORMAL  <- 1000   # cap normal reference cells too

tumour_cells <- colnames(seurat)[seurat$tumour_status == "tumour"]
normal_cells <- colnames(seurat)[seurat$tumour_status == "normal_reference"]

if (length(tumour_cells) > MAX_TUMOUR)
  tumour_cells <- sample(tumour_cells, MAX_TUMOUR)
if (length(normal_cells) > MAX_NORMAL)
  normal_cells <- sample(normal_cells, MAX_NORMAL)

keep_cells <- c(tumour_cells, normal_cells)
seurat     <- seurat[, keep_cells]
cat(sprintf("      Subsampled: %d tumour + %d normal = %d total cells\n",
            length(tumour_cells), length(normal_cells), ncol(seurat)))

# ── 2. Prepare InferCNV inputs ────────────────────────────────────────────────
cat("[2/4] Preparing InferCNV inputs...\n")

# Raw counts (InferCNV requires raw, not normalised)
raw_counts <- GetAssayData(seurat, assay = "RNA", layer = "counts")

# Cell annotations: tumour vs. normal_reference
ann_path <- file.path(OUT_DIR, "cell_annotations.txt")
annotations <- data.frame(
  row.names = colnames(seurat),
  group     = seurat$tumour_status
)
write.table(annotations, ann_path, sep = "\t", quote = FALSE, col.names = FALSE)
cat(sprintf("      Tumour cells: %d | Normal reference: %d\n",
            sum(seurat$tumour_status == "tumour"),
            sum(seurat$tumour_status == "normal_reference")))

# ── 3. Create InferCNV object ─────────────────────────────────────────────────
cat("[3/4] Creating InferCNV object...\n")
infercnv_obj <- CreateInfercnvObject(
  raw_counts_matrix = raw_counts,
  annotations_file  = ann_path,
  delim             = "\t",
  gene_order_file   = GENE_ORDER,
  ref_group_names   = "normal_reference"
)

# ── 4. Run InferCNV ───────────────────────────────────────────────────────────
cat("[4/4] Running InferCNV (this takes 20-60 min)...\n")
infercnv_obj <- infercnv::run(
  infercnv_obj,
  cutoff             = 0.1,      # min mean expression to include gene
  out_dir            = OUT_DIR,
  cluster_by_groups  = TRUE,     # cluster tumour + normal separately
  denoise            = TRUE,
  HMM                = TRUE,     # Hidden Markov Model for discrete CNV states
  HMM_type           = "i6",     # 6-state HMM: 0,0.5,1,1.5,2,3 copies
  num_threads        = 4,
  no_prelim_plot     = TRUE,     # skip intermediate plots to save time
  output_format      = "pdf"
)

# Save the object for later use in consensus step
saveRDS(infercnv_obj, file.path(OUT_DIR, "infercnv_obj.rds"))

# ── Extract per-cell CNV scores per chromosome arm ────────────────────────────
# The HMM states give discrete copy number per gene; average per arm
cat("[done] Extracting arm-level CNV scores...\n")

# Load chromosome arm annotations
arm_bed <- read.table(
  file.path(ROOT, "data", "reference", "chromosome_arms_hg38.bed"),
  col.names = c("chrom", "start", "end", "arm")
)

# Gene order file maps gene → chromosome position
gene_order <- read.table(GENE_ORDER, col.names = c("gene", "chrom", "start", "end"))
gene_order$arm <- paste0(
  gene_order$chrom,
  ifelse(gene_order$start < ave(gene_order$end[gene_order$chrom == gene_order$chrom],
                                 gene_order$chrom, FUN = median), "p", "q")
)

# Pull HMM CNV matrix (genes × cells)
cnv_mat <- infercnv_obj@expr.data

# Summarise to arm level (mean CNV signal per arm per cell)
arm_scores <- do.call(rbind, lapply(unique(gene_order$arm), function(a) {
  genes_in_arm <- intersect(gene_order$gene[gene_order$arm == a], rownames(cnv_mat))
  if (length(genes_in_arm) < 5) return(NULL)
  arm_mean <- colMeans(cnv_mat[genes_in_arm, , drop = FALSE])
  data.frame(arm = a, cell = names(arm_mean), cnv_score = arm_mean,
             row.names = NULL)
}))

write.csv(arm_scores,
          file.path(OUT_DIR, "infercnv_arm_scores.csv"),
          row.names = FALSE)
cat(sprintf("[done] Arm-level CNV scores saved → results/cnv_calls/infercnv/infercnv_arm_scores.csv\n"))
cat("[done] InferCNV complete. Main heatmap in results/cnv_calls/infercnv/\n")
cat("\nNext: Rscript pipeline/02_cnv/run_copykat.R\n")
