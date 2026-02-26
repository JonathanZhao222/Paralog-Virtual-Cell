# Phase 4 — Expression Validation
#
# For each high-confidence arm loss identified in Phase 3, confirm that genes on
# that arm show significantly lower expression in cells carrying the loss versus
# cells without it.
#
# Method:
#   - Split tumour cells into "arm-loss" vs. "no-loss" groups per arm.
#   - For each gene on the arm, run a Wilcoxon rank-sum test.
#   - Summarise as a volcano plot (log2 fold-change vs. -log10 p-value).
#   - Also produce a box plot of per-arm median expression for the top figures.
#
# Input:
#   data/processed/seurat_qc.rds
#   results/consensus/cell_cnv_status.csv
#   results/consensus/arm_loss_freq.csv
#   data/reference/gene_order_hg38.txt
#
# Output:
#   results/expression_validation/arm_de_results.csv
#   figures/04_expression/volcano_<arm>.pdf   (one per lost arm)
#   figures/04_expression/arm_expression_boxplot.pdf

library(Seurat)
library(ggplot2)
library(dplyr)

# ── Paths ──────────────────────────────────────────────────────────────────────
ROOT         <- getwd()
SEURAT_RDS   <- file.path(ROOT, "data", "processed", "seurat_qc.rds")
CNV_STATUS   <- file.path(ROOT, "results", "consensus", "cell_cnv_status.csv")
ARM_FREQ     <- file.path(ROOT, "results", "consensus", "arm_loss_freq.csv")
GENE_ORDER   <- file.path(ROOT, "data", "reference", "gene_order_hg38.txt")
ARM_BED      <- file.path(ROOT, "data", "reference", "chromosome_arms_hg38.bed")
OUT_DIR      <- file.path(ROOT, "results", "expression_validation")
FIG_DIR      <- file.path(ROOT, "figures", "04_expression")
dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)
dir.create(FIG_DIR, showWarnings = FALSE, recursive = TRUE)

# P-value and fold-change thresholds for volcano highlighting
PADJ_THRESH <- 0.05
LFC_THRESH  <- 0.5   # |log2 FC| threshold for "significant" labelling

cat("=== Phase 4: Expression Validation ===\n\n")

# ── 1. Load data ───────────────────────────────────────────────────────────────
cat("[1/4] Loading data...\n")
seurat     <- readRDS(SEURAT_RDS)
cell_status <- read.csv(CNV_STATUS)
arm_freq    <- read.csv(ARM_FREQ)
gene_order  <- read.table(GENE_ORDER, col.names = c("gene", "chrom", "start", "end"))
arm_bed     <- read.table(ARM_BED,    col.names = c("chrom", "start", "end", "arm"))

cat(sprintf("      %d cells in Seurat | %d high-confidence lost arms\n",
            ncol(seurat), nrow(arm_freq)))

if (nrow(arm_freq) == 0) {
  cat("      No high-confidence arm losses found. Exiting.\n")
  quit(save = "no", status = 0)
}

# ── 2. Assign chromosome arm to each gene ─────────────────────────────────────
cat("[2/4] Mapping genes to chromosome arms...\n")

assign_arm <- function(chrom, pos) {
  hits <- arm_bed[arm_bed$chrom == chrom & arm_bed$start <= pos & arm_bed$end > pos, ]
  if (nrow(hits) == 0) return(NA_character_)
  hits$arm[1]
}

gene_order$arm <- mapply(assign_arm, gene_order$chrom,
                         floor((gene_order$start + gene_order$end) / 2))
gene_order <- gene_order[!is.na(gene_order$arm), ]

# ── 3. Normalised expression matrix ───────────────────────────────────────────
cat("[3/4] Extracting normalised expression...\n")

# Use log-normalised counts (SCT or lognorm layer; fall back gracefully)
norm_layer <- tryCatch(
  GetAssayData(seurat, assay = "RNA", layer = "data"),
  error = function(e) GetAssayData(seurat, assay = "RNA", layer = "counts")
)

# ── 4. Per-arm differential expression ───────────────────────────────────────
cat("[4/4] Running Wilcoxon tests per arm...\n")

all_de <- list()

for (this_arm in arm_freq$lost_arm) {

  # Genes on this arm that are in the Seurat matrix
  arm_genes <- gene_order$gene[gene_order$arm == this_arm]
  arm_genes <- intersect(arm_genes, rownames(norm_layer))

  if (length(arm_genes) < 3) {
    cat(sprintf("      %s — skipped (only %d measured genes)\n",
                this_arm, length(arm_genes)))
    next
  }

  # Cells with arm loss (from the consensus status)
  loss_cells <- cell_status$cell[
    !is.na(cell_status$lost_arms) &
    grepl(this_arm, cell_status$lost_arms, fixed = TRUE)
  ]
  noloss_cells <- cell_status$cell[
    cell_status$is_tumour &
    !cell_status$cell %in% loss_cells
  ]

  # Keep only cells present in the Seurat object
  loss_cells   <- intersect(loss_cells,   colnames(norm_layer))
  noloss_cells <- intersect(noloss_cells, colnames(norm_layer))

  if (length(loss_cells) < 10 || length(noloss_cells) < 10) {
    cat(sprintf("      %s — skipped (loss: %d, no-loss: %d cells)\n",
                this_arm, length(loss_cells), length(noloss_cells)))
    next
  }

  cat(sprintf("      %s — %d genes, %d loss vs %d no-loss cells\n",
              this_arm, length(arm_genes), length(loss_cells), length(noloss_cells)))

  # Wilcoxon test per gene
  de_rows <- lapply(arm_genes, function(g) {
    x_loss   <- as.numeric(norm_layer[g, loss_cells])
    x_noloss <- as.numeric(norm_layer[g, noloss_cells])
    wt <- suppressWarnings(wilcox.test(x_loss, x_noloss, exact = FALSE))
    # log2 fold-change: mean(loss) - mean(no-loss) in log-space ≈ log2(ratio)
    lfc <- mean(x_loss) - mean(x_noloss)
    data.frame(
      gene    = g,
      arm     = this_arm,
      lfc     = lfc,
      p_value = wt$p.value,
      n_loss  = length(loss_cells),
      n_noloss = length(noloss_cells)
    )
  })
  de_arm <- do.call(rbind, de_rows)
  de_arm$padj <- p.adjust(de_arm$p_value, method = "BH")
  all_de[[this_arm]] <- de_arm

  # ── Volcano plot ──────────────────────────────────────────────────────────
  de_arm$sig <- with(de_arm,
    ifelse(padj < PADJ_THRESH & lfc < -LFC_THRESH, "down",
    ifelse(padj < PADJ_THRESH & lfc >  LFC_THRESH, "up", "ns")))

  # Label top 10 most significantly down-regulated
  label_genes <- de_arm %>%
    filter(sig == "down") %>%
    arrange(padj, lfc) %>%
    head(10)

  p_volc <- ggplot(de_arm, aes(x = lfc, y = -log10(padj), colour = sig)) +
    geom_point(alpha = 0.5, size = 1.2) +
    scale_colour_manual(values = c("down" = "#E74C3C", "up" = "#3498DB", "ns" = "grey60")) +
    geom_vline(xintercept = c(-LFC_THRESH, LFC_THRESH), linetype = "dashed", colour = "grey40") +
    geom_hline(yintercept = -log10(PADJ_THRESH), linetype = "dashed", colour = "grey40") +
    ggrepel::geom_text_repel(
      data = label_genes, aes(label = gene),
      size = 2.5, colour = "black", max.overlaps = 20
    ) +
    labs(
      title    = sprintf("Expression on %s: Loss vs No-Loss Tumour Cells", this_arm),
      subtitle = sprintf("%d arm-loss cells | %d no-loss cells | %d genes tested",
                         length(loss_cells), length(noloss_cells), length(arm_genes)),
      x        = "Log2 Fold-Change (loss − no-loss)",
      y        = "-log10(adj. p-value)",
      colour   = NULL
    ) +
    theme_classic()
  ggsave(file.path(FIG_DIR, sprintf("volcano_%s.pdf", gsub("chr", "", this_arm))),
         p_volc, width = 7, height = 5)
}

# ── Save combined DE table ────────────────────────────────────────────────────
if (length(all_de) > 0) {
  de_all <- do.call(rbind, all_de)
  write.csv(de_all, file.path(OUT_DIR, "arm_de_results.csv"), row.names = FALSE)
  cat(sprintf("\n[done] DE results: %d genes across %d arms\n",
              nrow(de_all), length(all_de)))
  cat(sprintf("       Saved → results/expression_validation/arm_de_results.csv\n"))

  # ── Summary box plot: median expression by arm-loss status ───────────────
  # Average expression of all genes on each arm, per cell
  summary_rows <- lapply(names(all_de), function(this_arm) {
    arm_genes <- gene_order$gene[gene_order$arm == this_arm]
    arm_genes <- intersect(arm_genes, rownames(norm_layer))
    if (length(arm_genes) == 0) return(NULL)

    loss_cells   <- cell_status$cell[
      !is.na(cell_status$lost_arms) &
      grepl(this_arm, cell_status$lost_arms, fixed = TRUE)
    ]
    noloss_cells <- cell_status$cell[
      cell_status$is_tumour & !cell_status$cell %in% loss_cells
    ]
    loss_cells   <- intersect(loss_cells,   colnames(norm_layer))
    noloss_cells <- intersect(noloss_cells, colnames(norm_layer))
    if (length(loss_cells) < 5 || length(noloss_cells) < 5) return(NULL)

    loss_med   <- median(colMeans(norm_layer[arm_genes, loss_cells,   drop = FALSE]))
    noloss_med <- median(colMeans(norm_layer[arm_genes, noloss_cells, drop = FALSE]))
    data.frame(arm = this_arm,
               group = c("loss", "no-loss"),
               median_expr = c(loss_med, noloss_med))
  })
  summary_df <- do.call(rbind, Filter(Negate(is.null), summary_rows))

  if (nrow(summary_df) > 0) {
    p_box <- ggplot(summary_df, aes(x = arm, y = median_expr, fill = group)) +
      geom_col(position = "dodge") +
      scale_fill_manual(values = c("loss" = "#E74C3C", "no-loss" = "#3498DB")) +
      labs(
        title = "Median Arm Expression: Loss vs No-Loss Tumour Cells",
        x     = "Chromosome Arm",
        y     = "Median normalised expression (per cell)",
        fill  = NULL
      ) +
      theme_classic() +
      theme(axis.text.x = element_text(angle = 45, hjust = 1))
    ggsave(file.path(FIG_DIR, "arm_expression_boxplot.pdf"), p_box, width = 10, height = 5)
    cat("       Saved → figures/04_expression/arm_expression_boxplot.pdf\n")
  }
} else {
  cat("\n[done] No arms had enough cells for testing.\n")
}

cat("\nNext: Rscript pipeline/05_paralog_analysis/query_ensembl_paralogs.R\n")
