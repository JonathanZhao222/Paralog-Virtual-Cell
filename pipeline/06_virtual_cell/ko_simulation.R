# Phase 6b — Paralog Dependency Scoring (Expression-Based)
#
# Rather than running a full GRN-based virtual KO (scTenifoldKnk), we compute a
# combined expression-based dependency score from the Phase 5b results.
# This is computationally efficient and directly interpretable.
#
# KO Dependency Score = LFC_paralog × (−Spearman_r)
#
# Intuition:
#   LFC_paralog > 0 : paralog P is upregulated in cells that have lost G's arm
#   Spearman_r < 0  : P expression anti-correlates with G in loss cells
#                     (when G is low, P is high → compensation / buffering)
#   High combined score → P is both induced AND anti-correlated with G
#                       → strong evidence of paralog dependency
#
# Additionally computes a "virtual KO" approximation:
#   In arm-loss cells, stratify by G expression quartile.
#   Compare P expression in lowest-G quartile (G virtually absent)
#   vs highest-G quartile (G present). A large difference = compensation.
#
# Input:
#   results/paralog_pairs/paralog_upregulation.csv
#   results/paralog_pairs/paralog_correlation.csv
#   results/consensus/cell_cnv_status.csv
#   data/processed/seurat_qc.rds
#
# Output:
#   results/ko_scores/ko_scores.csv
#   results/ko_scores/ko_arm_summary.csv
#   figures/06_virtual_cell/ko_score_dotplot.pdf
#   figures/06_virtual_cell/top_ko_boxplots.pdf

library(Seurat)
library(ggplot2)
library(dplyr)

# ── Paths ──────────────────────────────────────────────────────────────────────
ROOT        <- getwd()
UPREG       <- file.path(ROOT, "results", "paralog_pairs", "paralog_upregulation.csv")
CORR        <- file.path(ROOT, "results", "paralog_pairs", "paralog_correlation.csv")
CNV_STATUS  <- file.path(ROOT, "results", "consensus", "cell_cnv_status.csv")
SEURAT_RDS  <- file.path(ROOT, "data", "processed", "seurat_qc.rds")
OUT_DIR     <- file.path(ROOT, "results", "ko_scores")
FIG_DIR     <- file.path(ROOT, "figures", "06_virtual_cell")
dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)
dir.create(FIG_DIR, showWarnings = FALSE, recursive = TRUE)

cat("=== Phase 6b: Paralog Dependency Scoring ===\n\n")

# ── 1. Load Phase 5b results ───────────────────────────────────────────────────
cat("[1/4] Loading Phase 5b results...\n")
upreg       <- read.csv(UPREG)
corr        <- read.csv(CORR)
cell_status <- read.csv(CNV_STATUS)
seurat      <- readRDS(SEURAT_RDS)

norm_mat <- tryCatch(
  GetAssayData(seurat, assay = "RNA", layer = "data"),
  error = function(e) GetAssayData(seurat, assay = "RNA", layer = "counts")
)

cat(sprintf("      Upregulation table: %d pairs\n", nrow(upreg)))
cat(sprintf("      Correlation table:  %d pairs\n", nrow(corr)))

# ── 2. Compute combined KO dependency score ───────────────────────────────────
cat("[2/4] Computing KO dependency scores...\n")

# Join upregulation and correlation tables
ko_df <- upreg %>%
  inner_join(
    corr %>% select(gene, paralog_gene, query_arm, spearman_r, padj) %>%
      rename(spearman_padj = padj),
    by = c("gene", "paralog_gene", "query_arm")
  ) %>%
  mutate(
    # KO dependency score: high when P is upregulated AND anti-correlated with G
    # Clip spearman_r to avoid sign issues from noise near 0
    ko_score = lfc_paralog * pmax(-spearman_r, 0)
  ) %>%
  arrange(desc(ko_score))

cat(sprintf("      Scored pairs: %d\n", nrow(ko_df)))

# ── 3. Virtual KO approximation: G-low vs G-high quartile comparison ──────────
cat("[3/4] Computing virtual KO approximation (G-low vs G-high quartile)...\n")

top_pairs <- ko_df %>%
  filter(padj < 0.05, lfc_paralog > 0.3, ko_score > 0) %>%
  head(50)

vko_rows <- list()
for (idx in seq_len(nrow(top_pairs))) {
  G   <- top_pairs$gene[idx]
  P   <- top_pairs$paralog_gene[idx]
  arm <- top_pairs$query_arm[idx]

  if (!G %in% rownames(norm_mat) || !P %in% rownames(norm_mat)) next

  loss_cells <- cell_status$cell[
    !is.na(cell_status$lost_arms) &
    grepl(arm, cell_status$lost_arms, fixed = TRUE)
  ]
  loss_cells <- intersect(loss_cells, colnames(norm_mat))
  if (length(loss_cells) < 20) next

  g_expr <- as.numeric(norm_mat[G, loss_cells])
  p_expr <- as.numeric(norm_mat[P, loss_cells])

  # Stratify by G expression quartile
  q25 <- quantile(g_expr, 0.25)
  q75 <- quantile(g_expr, 0.75)
  low_g  <- p_expr[g_expr <= q25]
  high_g <- p_expr[g_expr >= q75]

  if (length(low_g) < 5 || length(high_g) < 5) next

  wt  <- suppressWarnings(wilcox.test(low_g, high_g, exact = FALSE))
  vko_rows[[idx]] <- data.frame(
    gene              = G,
    paralog_gene      = P,
    query_arm         = arm,
    mean_P_g_low      = mean(low_g),
    mean_P_g_high     = mean(high_g),
    vko_lfc           = mean(low_g) - mean(high_g),  # P higher when G is low?
    vko_pvalue        = wt$p.value,
    n_low_cells       = length(low_g),
    n_high_cells      = length(high_g)
  )
}
vko_df <- do.call(rbind, Filter(Negate(is.null), vko_rows))

if (!is.null(vko_df) && nrow(vko_df) > 0) {
  vko_df$vko_padj <- p.adjust(vko_df$vko_pvalue, method = "BH")
  ko_df <- ko_df %>%
    left_join(vko_df %>% select(gene, paralog_gene, query_arm,
                                 vko_lfc, vko_padj),
              by = c("gene", "paralog_gene", "query_arm"))
}

write.csv(ko_df, file.path(OUT_DIR, "ko_scores.csv"), row.names = FALSE)
cat(sprintf("      Saved → results/ko_scores/ko_scores.csv\n"))

# Top result
top1 <- ko_df[1, ]
cat(sprintf("      Top dependency: %s → %s  score=%.3f  arm=%s\n",
            top1$gene, top1$paralog_gene, top1$ko_score, top1$query_arm))

# ── 4. Figures ─────────────────────────────────────────────────────────────────
cat("[4/4] Generating figures...\n")

# ── Dot plot: top 30 pairs ──────────────────────────────────────────────────
top30 <- ko_df %>%
  filter(padj < 0.05, lfc_paralog > 0.3) %>%
  head(30)

if (nrow(top30) > 0) {
  top30$label <- paste(top30$gene, "\u2192", top30$paralog_gene)
  top30$label <- factor(top30$label, levels = rev(top30$label))
  top30$neg_log_padj <- pmin(-log10(top30$padj + 1e-300), 15)

  p_dot <- ggplot(top30,
                  aes(x = ko_score, y = label,
                      colour = ko_score, size = neg_log_padj)) +
    geom_point() +
    scale_colour_gradient(low = "#F39C12", high = "#C0392B") +
    scale_size_continuous(range = c(2, 8), name = "-log10(padj)") +
    labs(
      title    = "Paralog Dependency Scores",
      subtitle = "Score = paralog upregulation × anti-correlation with lost gene",
      x        = "KO Dependency Score",
      y        = "Lost gene \u2192 Paralog",
      colour   = "Score"
    ) +
    theme_classic() +
    theme(axis.text.y = element_text(size = 8))
  ggsave(file.path(FIG_DIR, "ko_score_dotplot.pdf"), p_dot,
         width = 9, height = max(6, nrow(top30) * 0.35))
  cat("      Saved → figures/06_virtual_cell/ko_score_dotplot.pdf\n")
}

# ── Virtual KO box plots: top 6 pairs ─────────────────────────────────────
if (!is.null(vko_df) && nrow(vko_df) > 0) {
  top6 <- vko_df %>%
    filter(vko_padj < 0.05, vko_lfc > 0) %>%
    arrange(vko_padj) %>%
    head(6)

  if (nrow(top6) > 0) {
    box_plots <- lapply(seq_len(nrow(top6)), function(i) {
      G   <- top6$gene[i]
      P   <- top6$paralog_gene[i]
      arm <- top6$query_arm[i]

      loss_cells <- cell_status$cell[
        !is.na(cell_status$lost_arms) &
        grepl(arm, cell_status$lost_arms, fixed = TRUE)
      ]
      loss_cells <- intersect(loss_cells, colnames(norm_mat))
      if (!G %in% rownames(norm_mat) || !P %in% rownames(norm_mat)) return(NULL)

      g_expr <- as.numeric(norm_mat[G, loss_cells])
      p_expr <- as.numeric(norm_mat[P, loss_cells])
      q25 <- quantile(g_expr, 0.25)
      q75 <- quantile(g_expr, 0.75)

      df_box <- data.frame(
        P_expr = c(p_expr[g_expr <= q25], p_expr[g_expr >= q75]),
        G_level = c(rep("G low (virtual KO)", sum(g_expr <= q25)),
                    rep("G high (control)", sum(g_expr >= q75)))
      )
      df_box$G_level <- factor(df_box$G_level,
                                levels = c("G low (virtual KO)", "G high (control)"))

      ggplot(df_box, aes(x = G_level, y = P_expr, fill = G_level)) +
        geom_boxplot(outlier.size = 0.5, alpha = 0.8) +
        scale_fill_manual(values = c("G low (virtual KO)" = "#E74C3C",
                                     "G high (control)"   = "#3498DB")) +
        labs(
          title    = sprintf("%s \u2192 %s", G, P),
          subtitle = sprintf("arm: %s | vKO padj=%.3g", arm, top6$vko_padj[i]),
          x = NULL, y = sprintf("%s expression", P)
        ) +
        theme_classic(base_size = 9) +
        theme(legend.position = "none")
    })
    box_plots <- Filter(Negate(is.null), box_plots)

    if (length(box_plots) > 0) {
      library(patchwork)
      n_col <- min(3, length(box_plots))
      p_box <- wrap_plots(box_plots, ncol = n_col)
      ggsave(file.path(FIG_DIR, "top_ko_boxplots.pdf"), p_box,
             width = n_col * 3.5, height = ceiling(length(box_plots) / n_col) * 3.5)
      cat("      Saved → figures/06_virtual_cell/top_ko_boxplots.pdf\n")
    }
  }
}

# ── Arm-level summary ─────────────────────────────────────────────────────────
arm_summary <- ko_df %>%
  filter(padj < 0.05, lfc_paralog > 0.3) %>%
  group_by(query_arm) %>%
  summarise(
    n_sig_pairs    = n(),
    mean_ko_score  = mean(ko_score, na.rm = TRUE),
    top_gene       = gene[which.max(ko_score)],
    top_paralog    = paralog_gene[which.max(ko_score)],
    top_score      = max(ko_score, na.rm = TRUE),
    .groups        = "drop"
  ) %>%
  arrange(desc(mean_ko_score))

print(arm_summary)
write.csv(arm_summary, file.path(OUT_DIR, "ko_arm_summary.csv"), row.names = FALSE)
cat("      Saved → results/ko_scores/ko_arm_summary.csv\n")

cat("\n[done] Phase 6b complete.\n")
cat("\n=== PIPELINE COMPLETE ===\n")
cat("Key output figures:\n")
cat("  figures/02_cnv/copykat_umap.pdf\n")
cat("  figures/03_consensus/consensus_umap.pdf\n")
cat("  figures/03_consensus/arm_loss_barplot.pdf\n")
cat("  figures/04_expression/arm_expression_boxplot.pdf\n")
cat("  figures/04_expression/volcano_*.pdf\n")
cat("  figures/05_paralog/paralog_volcano.pdf\n")
cat("  figures/05_paralog/top_pairs_scatter.pdf\n")
cat("  figures/06_virtual_cell/ko_score_dotplot.pdf\n")
cat("  figures/06_virtual_cell/top_ko_boxplots.pdf\n")
