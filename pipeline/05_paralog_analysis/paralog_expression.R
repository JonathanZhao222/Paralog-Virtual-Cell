# Phase 5b — Paralog Expression Analysis
#
# For each cross-arm paralog pair (gene G on lost arm, paralog P on intact arm):
#
#   1. Upregulation test: is P more highly expressed in arm-loss cells than
#      no-loss tumour cells? (compensation signal)
#   2. Correlation test: in arm-loss cells, does P expression anti-correlate
#      with G expression? (when G is low, P is high → buffering)
#   3. Summarise as a scatter plot (G vs P expression coloured by arm-loss status)
#      and a summary volcano of paralog upregulation.
#
# Input:
#   data/processed/seurat_qc.rds
#   results/consensus/cell_cnv_status.csv
#   results/paralog_pairs/cross_arm_paralogs.csv
#
# Output:
#   results/paralog_pairs/paralog_upregulation.csv   — per-pair Wilcoxon results
#   results/paralog_pairs/paralog_correlation.csv    — per-pair Spearman r
#   figures/05_paralog/paralog_volcano.pdf
#   figures/05_paralog/top_pairs_scatter.pdf

library(Seurat)
library(ggplot2)
library(dplyr)

# ── Paths ──────────────────────────────────────────────────────────────────────
ROOT       <- getwd()
SEURAT_RDS <- file.path(ROOT, "data", "processed", "seurat_qc.rds")
CNV_STATUS <- file.path(ROOT, "results", "consensus", "cell_cnv_status.csv")
PARALOGS   <- file.path(ROOT, "results", "paralog_pairs", "cross_arm_paralogs.csv")
OUT_DIR    <- file.path(ROOT, "results", "paralog_pairs")
FIG_DIR    <- file.path(ROOT, "figures", "05_paralog")
dir.create(FIG_DIR, showWarnings = FALSE, recursive = TRUE)

PADJ_THRESH <- 0.05
LFC_THRESH  <- 0.3
MIN_CELLS   <- 20   # minimum cells per group for testing

cat("=== Phase 5b: Paralog Expression Analysis ===\n\n")

# ── 1. Load data ───────────────────────────────────────────────────────────────
cat("[1/3] Loading data...\n")
seurat      <- readRDS(SEURAT_RDS)
cell_status <- read.csv(CNV_STATUS)
paralogs    <- read.csv(PARALOGS)

cat(sprintf("      %d cross-arm paralog pairs\n", nrow(paralogs)))

norm_mat <- tryCatch(
  GetAssayData(seurat, assay = "RNA", layer = "data"),
  error = function(e) GetAssayData(seurat, assay = "RNA", layer = "counts")
)

# ── 2. Paralog upregulation test ───────────────────────────────────────────────
cat("[2/3] Testing paralog upregulation in arm-loss cells...\n")

upregulation_rows <- list()
correlation_rows  <- list()

# Work per query_arm so we use the right loss/no-loss cell sets
for (this_arm in unique(paralogs$query_arm)) {

  loss_cells   <- cell_status$cell[
    !is.na(cell_status$lost_arms) &
    grepl(this_arm, cell_status$lost_arms, fixed = TRUE)
  ]
  noloss_cells <- cell_status$cell[
    cell_status$is_tumour & !cell_status$cell %in% loss_cells
  ]
  loss_cells   <- intersect(loss_cells,   colnames(norm_mat))
  noloss_cells <- intersect(noloss_cells, colnames(norm_mat))

  if (length(loss_cells) < MIN_CELLS || length(noloss_cells) < MIN_CELLS) next

  arm_pairs <- paralogs[paralogs$query_arm == this_arm, ]
  # Keep only pairs where BOTH genes are measured
  arm_pairs <- arm_pairs[
    arm_pairs$gene        %in% rownames(norm_mat) &
    arm_pairs$paralog_gene %in% rownames(norm_mat),
  ]

  if (nrow(arm_pairs) == 0) next

  cat(sprintf("      %s: %d pairs, %d loss vs %d no-loss cells\n",
              this_arm, nrow(arm_pairs), length(loss_cells), length(noloss_cells)))

  for (idx in seq_len(nrow(arm_pairs))) {
    G <- arm_pairs$gene[idx]
    P <- arm_pairs$paralog_gene[idx]

    p_loss   <- as.numeric(norm_mat[P, loss_cells])
    p_noloss <- as.numeric(norm_mat[P, noloss_cells])
    g_loss   <- as.numeric(norm_mat[G, loss_cells])

    wt   <- suppressWarnings(wilcox.test(p_loss, p_noloss, exact = FALSE))
    lfc_P <- mean(p_loss) - mean(p_noloss)

    # Spearman correlation between G and P expression in arm-loss cells
    sp <- suppressWarnings(cor.test(g_loss, p_loss, method = "spearman", exact = FALSE))

    upregulation_rows[[length(upregulation_rows) + 1]] <- data.frame(
      query_arm    = this_arm,
      gene         = G,
      paralog_gene = P,
      lfc_paralog  = lfc_P,
      p_value      = wt$p.value,
      n_loss       = length(loss_cells),
      n_noloss     = length(noloss_cells),
      mean_G_loss  = mean(g_loss),
      mean_P_loss  = mean(p_loss),
      mean_P_noloss = mean(p_noloss)
    )

    correlation_rows[[length(correlation_rows) + 1]] <- data.frame(
      query_arm    = this_arm,
      gene         = G,
      paralog_gene = P,
      spearman_r   = sp$estimate,
      spearman_p   = sp$p.value,
      n_cells      = length(loss_cells)
    )
  }
}

# Combine and adjust p-values
upreg_df <- do.call(rbind, upregulation_rows)
corr_df  <- do.call(rbind, correlation_rows)

if (!is.null(upreg_df) && nrow(upreg_df) > 0) {
  upreg_df$padj <- p.adjust(upreg_df$p_value, method = "BH")
  write.csv(upreg_df, file.path(OUT_DIR, "paralog_upregulation.csv"), row.names = FALSE)
  cat(sprintf("\n      Upregulation table: %d pairs tested\n", nrow(upreg_df)))
  cat(sprintf("      Significant (padj < %.2f, LFC > %.1f): %d pairs\n",
              PADJ_THRESH, LFC_THRESH,
              sum(upreg_df$padj < PADJ_THRESH & upreg_df$lfc_paralog > LFC_THRESH,
                  na.rm = TRUE)))
  cat("      Saved → results/paralog_pairs/paralog_upregulation.csv\n")
}

if (!is.null(corr_df) && nrow(corr_df) > 0) {
  corr_df$padj <- p.adjust(corr_df$spearman_p, method = "BH")
  write.csv(corr_df, file.path(OUT_DIR, "paralog_correlation.csv"), row.names = FALSE)
  cat("      Saved → results/paralog_pairs/paralog_correlation.csv\n")
}

# ── 3. Figures ─────────────────────────────────────────────────────────────────
cat("[3/3] Generating figures...\n")

# ── Volcano: paralog upregulation ────────────────────────────────────────────
if (!is.null(upreg_df) && nrow(upreg_df) > 0) {
  upreg_df$sig <- with(upreg_df,
    ifelse(padj < PADJ_THRESH & lfc_paralog >  LFC_THRESH, "upregulated",
    ifelse(padj < PADJ_THRESH & lfc_paralog < -LFC_THRESH, "downregulated", "ns")))

  label_pairs <- upreg_df %>%
    filter(sig == "upregulated") %>%
    arrange(padj, desc(lfc_paralog)) %>%
    head(15) %>%
    mutate(label = paste(gene, "→", paralog_gene))

  p_volc <- ggplot(upreg_df,
                   aes(x = lfc_paralog, y = -log10(padj), colour = sig)) +
    geom_point(alpha = 0.5, size = 1.2) +
    scale_colour_manual(values = c("upregulated"   = "#27AE60",
                                   "downregulated" = "#E74C3C",
                                   "ns"            = "grey70")) +
    geom_vline(xintercept = c(-LFC_THRESH, LFC_THRESH),
               linetype = "dashed", colour = "grey40") +
    geom_hline(yintercept = -log10(PADJ_THRESH),
               linetype = "dashed", colour = "grey40") +
    ggrepel::geom_text_repel(
      data = label_pairs,
      aes(label = label), size = 2.5, colour = "black", max.overlaps = 20
    ) +
    labs(
      title    = "Paralog Upregulation in Arm-Loss Cells",
      subtitle = "Green = paralog significantly UP in cells that have lost the query gene's arm",
      x        = "Log2 FC (arm-loss vs no-loss) for PARALOG",
      y        = "-log10(adj. p-value)",
      colour   = NULL
    ) +
    theme_classic()
  ggsave(file.path(FIG_DIR, "paralog_volcano.pdf"), p_volc, width = 8, height = 6)
  cat("      Saved → figures/05_paralog/paralog_volcano.pdf\n")
}

# ── Top pair scatter plots ────────────────────────────────────────────────────
# For the top 6 compensating paralog pairs, plot G vs P in loss / no-loss cells
if (!is.null(upreg_df) && nrow(upreg_df) > 0) {
  top_pairs <- upreg_df %>%
    filter(padj < PADJ_THRESH, lfc_paralog > LFC_THRESH) %>%
    arrange(padj, desc(lfc_paralog)) %>%
    head(6)

  if (nrow(top_pairs) > 0) {
    scatter_plots <- lapply(seq_len(nrow(top_pairs)), function(i) {
      G   <- top_pairs$gene[i]
      P   <- top_pairs$paralog_gene[i]
      arm <- top_pairs$query_arm[i]

      loss_cells <- cell_status$cell[
        !is.na(cell_status$lost_arms) &
        grepl(arm, cell_status$lost_arms, fixed = TRUE)
      ]
      noloss_cells <- cell_status$cell[
        cell_status$is_tumour & !cell_status$cell %in% loss_cells
      ]
      all_test_cells <- intersect(c(loss_cells, noloss_cells), colnames(norm_mat))
      all_test_cells <- all_test_cells[
        all_test_cells %in% colnames(norm_mat) &
        G %in% rownames(norm_mat) & P %in% rownames(norm_mat)
      ]
      if (length(all_test_cells) < 10) return(NULL)

      df_sc <- data.frame(
        G     = as.numeric(norm_mat[G, all_test_cells]),
        P     = as.numeric(norm_mat[P, all_test_cells]),
        group = ifelse(all_test_cells %in% loss_cells, "arm-loss", "no-loss")
      )

      ggplot(df_sc, aes(x = G, y = P, colour = group)) +
        geom_point(alpha = 0.3, size = 0.8) +
        geom_smooth(method = "lm", se = FALSE, linewidth = 0.8) +
        scale_colour_manual(values = c("arm-loss" = "#E74C3C", "no-loss" = "#3498DB")) +
        labs(
          title  = sprintf("%s (lost) → %s (paralog)", G, P),
          subtitle = sprintf("Arm: %s  |  padj = %.3g  |  LFC = %.2f",
                              arm, top_pairs$padj[i], top_pairs$lfc_paralog[i]),
          x      = sprintf("%s expression (log-norm)", G),
          y      = sprintf("%s expression (log-norm)", P),
          colour = NULL
        ) +
        theme_classic(base_size = 9)
    })
    scatter_plots <- Filter(Negate(is.null), scatter_plots)

    if (length(scatter_plots) > 0) {
      library(patchwork)
      n_col <- min(3, length(scatter_plots))
      p_scatter <- wrap_plots(scatter_plots, ncol = n_col)
      ggsave(file.path(FIG_DIR, "top_pairs_scatter.pdf"), p_scatter,
             width = n_col * 4, height = ceiling(length(scatter_plots) / n_col) * 3.5)
      cat("      Saved → figures/05_paralog/top_pairs_scatter.pdf\n")
    }
  }
}

cat("\n[done] Phase 5b complete.\n")
cat("\nNext: Rscript pipeline/06_virtual_cell/sctenifold_network.R\n")
