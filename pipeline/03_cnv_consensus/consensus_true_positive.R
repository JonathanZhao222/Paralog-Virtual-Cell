# Phase 3 — CNV Consensus: True-Positive Arm Losses
#
# Strategy:
#   1. InferCNV assigns HMM CNV states to tumour SUBCLUSTERS (not individual cells).
#      We map cells → their subcluster → arm-level loss calls.
#   2. CopyKAT assigns aneuploid / diploid to individual cells.
#   3. A tumour cell is called a "consensus aneuploid" if:
#        - CopyKAT calls it aneuploid, AND
#        - Its InferCNV subcluster has ≥ 1 arm with a loss state.
#   4. An arm loss is "high-confidence" when ≥ 20 % of consensus aneuploid cells
#      carry it (tunable via LOSS_FREQ_THRESH).
#
# Input:
#   results/cnv_calls/infercnv/  (InferCNV HMM outputs)
#   results/cnv_calls/copykat/copykat_aneuploid_calls.csv
#   data/reference/chromosome_arms_hg38.bed
#   data/processed/seurat_qc.rds
#
# Output:
#   results/consensus/cell_cnv_status.csv   — per cell: consensus call + lost arms
#   results/consensus/arm_loss_freq.csv     — per arm: fraction of aneuploidy cells
#   figures/03_consensus/consensus_umap.pdf
#   figures/03_consensus/arm_loss_barplot.pdf

library(Seurat)
library(ggplot2)
library(dplyr)

# ── Paths ──────────────────────────────────────────────────────────────────────
ROOT        <- getwd()
INFERCNV_DIR <- file.path(ROOT, "results", "cnv_calls", "infercnv")
COPYKAT_DIR  <- file.path(ROOT, "results", "cnv_calls", "copykat")
ARM_BED     <- file.path(ROOT, "data", "reference", "chromosome_arms_hg38.bed")
SEURAT_RDS  <- file.path(ROOT, "data", "processed", "seurat_qc.rds")
OUT_DIR     <- file.path(ROOT, "results", "consensus")
FIG_DIR     <- file.path(ROOT, "figures", "03_consensus")
dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)
dir.create(FIG_DIR, showWarnings = FALSE, recursive = TRUE)

# Tuning parameters
LOSS_STATE_MAX   <- 3    # HMM states 1-3 = losses (0, 0.5, 1 copy)
GAIN_STATE_MIN   <- 6    # HMM state 6   = gain (3+ copies)
LOSS_FREQ_THRESH <- 0.20 # arm must be lost in ≥ 20 % of aneuploid cells

cat("=== Phase 3: CNV Consensus ===\n\n")

# ── 1. Load data ───────────────────────────────────────────────────────────────
cat("[1/6] Loading input data...\n")

# Chromosome arm boundaries
arm_bed <- read.table(ARM_BED, col.names = c("chrom", "start", "end", "arm"))

# InferCNV: cell → subcluster mapping
cell_groupings_path <- file.path(
  INFERCNV_DIR,
  "17_HMM_predHMMi6.leiden.hmm_mode-subclusters.cell_groupings"
)
cell_groups <- read.table(cell_groupings_path, header = TRUE, sep = "\t",
                          col.names = c("cell_group_name", "cell"))
cat(sprintf("      InferCNV: %d cells in %d subclusters\n",
            nrow(cell_groups), length(unique(cell_groups$cell_group_name))))

# InferCNV: subcluster → CNV regions (arm-level HMM states)
regions_path <- file.path(
  INFERCNV_DIR,
  "17_HMM_predHMMi6.leiden.hmm_mode-subclusters.pred_cnv_regions.dat"
)
cnv_regions <- read.table(regions_path, header = TRUE, sep = "\t")
cat(sprintf("      InferCNV: %d CNV region calls\n", nrow(cnv_regions)))

# CopyKAT: per-cell aneuploid/diploid calls
copykat_path <- file.path(COPYKAT_DIR, "copykat_aneuploid_calls.csv")
copykat <- read.csv(copykat_path)
cat(sprintf("      CopyKAT: %d cells (%d aneuploid, %d diploid)\n",
            nrow(copykat),
            sum(copykat$copykat_call == "aneuploid", na.rm = TRUE),
            sum(copykat$copykat_call == "diploid",   na.rm = TRUE)))

# Seurat object (for UMAP coordinates)
seurat <- readRDS(SEURAT_RDS)

# ── 2. Map InferCNV regions → chromosome arms ─────────────────────────────────
cat("[2/6] Mapping InferCNV regions to chromosome arms...\n")

assign_arm <- function(chrom, pos_start, pos_end) {
  # Use midpoint of CNV region
  mid <- (pos_start + pos_end) / 2
  hits <- arm_bed[arm_bed$chrom == chrom & arm_bed$start <= mid & arm_bed$end > mid, ]
  if (nrow(hits) == 0) return(NA_character_)
  hits$arm[1]
}

cnv_regions$arm <- mapply(
  assign_arm,
  cnv_regions$chr,
  cnv_regions$start,
  cnv_regions$end
)
cnv_regions <- cnv_regions[!is.na(cnv_regions$arm), ]
cat(sprintf("      %d regions mapped to %d unique arms\n",
            nrow(cnv_regions), length(unique(cnv_regions$arm))))

# ── 3. Arm-level loss per subcluster ─────────────────────────────────────────
cat("[3/6] Computing arm-level losses per InferCNV subcluster...\n")

# For each (subcluster, arm): classify as loss / normal / gain
# A subcluster-arm gets "loss" if ANY region in that arm has state ≤ LOSS_STATE_MAX
subcluster_arm_state <- cnv_regions %>%
  group_by(cell_group_name, arm) %>%
  summarise(
    min_state = min(state, na.rm = TRUE),
    max_state = max(state, na.rm = TRUE),
    .groups   = "drop"
  ) %>%
  mutate(
    arm_call = case_when(
      min_state <= LOSS_STATE_MAX ~ "loss",
      max_state >= GAIN_STATE_MIN ~ "gain",
      TRUE                        ~ "neutral"
    )
  )

# List of arms with loss per subcluster
subcluster_lost_arms <- subcluster_arm_state %>%
  filter(arm_call == "loss") %>%
  group_by(cell_group_name) %>%
  summarise(lost_arms = paste(arm, collapse = ";"), .groups = "drop")

# ── 4. Map cell → subcluster → lost arms ─────────────────────────────────────
cat("[4/6] Mapping individual cells to arm-level calls...\n")

cell_arm <- cell_groups %>%
  left_join(subcluster_lost_arms, by = "cell_group_name") %>%
  mutate(
    infercnv_has_loss = !is.na(lost_arms) & nchar(lost_arms) > 0,
    lost_arms         = ifelse(is.na(lost_arms), "", lost_arms)
  )

# Only tumour cells can be aneuploid
tumour_cells_all <- colnames(seurat)[seurat$tumour_status == "tumour"]
cell_arm$is_tumour <- cell_arm$cell %in% tumour_cells_all

# ── 5. Consensus call ─────────────────────────────────────────────────────────
cat("[5/6] Computing consensus aneuploid calls...\n")

# Join CopyKAT predictions
cell_status <- cell_arm %>%
  left_join(copykat %>% rename(cell = cell, copykat_call = copykat_call),
            by = "cell") %>%
  mutate(
    copykat_aneuploid = (!is.na(copykat_call) & copykat_call == "aneuploid"),
    consensus_aneuploid = is_tumour & infercnv_has_loss & copykat_aneuploid
  )

n_consensus <- sum(cell_status$consensus_aneuploid, na.rm = TRUE)
n_tumour    <- sum(cell_status$is_tumour, na.rm = TRUE)
cat(sprintf("      Tumour cells:      %d\n", n_tumour))
cat(sprintf("      Consensus aneuploid: %d (%.1f %% of tumour)\n",
            n_consensus, 100 * n_consensus / max(n_tumour, 1)))

write.csv(cell_status, file.path(OUT_DIR, "cell_cnv_status.csv"), row.names = FALSE)
cat(sprintf("      Saved → results/consensus/cell_cnv_status.csv\n"))

# ── 6. Arm-loss frequency across consensus aneuploid cells ────────────────────
cat("[6/6] Computing arm-loss frequencies and generating figures...\n")

aneuploid_cells <- cell_status %>% filter(consensus_aneuploid)

# Expand lost_arms column (semicolon-separated)
arm_loss_long <- aneuploid_cells %>%
  filter(nchar(lost_arms) > 0) %>%
  mutate(arm_list = strsplit(lost_arms, ";")) %>%
  tidyr::unnest(arm_list) %>%
  rename(lost_arm = arm_list)

arm_loss_freq <- arm_loss_long %>%
  count(lost_arm, name = "n_cells_with_loss") %>%
  mutate(
    total_aneuploid   = n_consensus,
    freq              = n_cells_with_loss / total_aneuploid
  ) %>%
  arrange(desc(freq)) %>%
  filter(freq >= LOSS_FREQ_THRESH)

cat(sprintf("      %d arms pass ≥%.0f%% frequency threshold\n",
            nrow(arm_loss_freq), 100 * LOSS_FREQ_THRESH))
write.csv(arm_loss_freq, file.path(OUT_DIR, "arm_loss_freq.csv"), row.names = FALSE)
cat(sprintf("      Saved → results/consensus/arm_loss_freq.csv\n"))

# ── Figure 1: Arm-loss bar chart ─────────────────────────────────────────────
if (nrow(arm_loss_freq) > 0) {
  arm_loss_freq$lost_arm <- factor(arm_loss_freq$lost_arm,
                                   levels = arm_loss_freq$lost_arm)
  p_bar <- ggplot(arm_loss_freq, aes(x = lost_arm, y = freq * 100)) +
    geom_col(fill = "#E74C3C") +
    geom_hline(yintercept = LOSS_FREQ_THRESH * 100, linetype = "dashed", colour = "grey40") +
    labs(
      title    = "Arm-level Loss Frequency (Consensus)",
      subtitle = sprintf("Among %d consensus-aneuploid tumour cells", n_consensus),
      x        = "Chromosome Arm",
      y        = "% Cells with Arm Loss"
    ) +
    theme_classic() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
  ggsave(file.path(FIG_DIR, "arm_loss_barplot.pdf"), p_bar, width = 10, height = 5)
  cat(sprintf("      Saved → figures/03_consensus/arm_loss_barplot.pdf\n"))
} else {
  cat("      WARNING: No arms passed the frequency threshold — barplot skipped.\n")
}

# ── Figure 2: UMAP coloured by consensus call ─────────────────────────────────
seurat$consensus_call <- "normal_reference"
# Cells in the subsampled InferCNV run — mark by available data
seurat$consensus_call[colnames(seurat) %in%
  cell_status$cell[cell_status$consensus_aneuploid]] <- "aneuploid"
seurat$consensus_call[colnames(seurat) %in%
  cell_status$cell[cell_status$is_tumour & !cell_status$consensus_aneuploid]] <- "tumour_unclassified"

p_umap <- DimPlot(seurat, reduction = "umap", group.by = "consensus_call",
                  cols = c(
                    "aneuploid"          = "#E74C3C",
                    "tumour_unclassified" = "#F39C12",
                    "normal_reference"   = "#3498DB"
                  )) +
  ggtitle("Consensus CNV: Aneuploid Tumour Cells") +
  theme_classic()
ggsave(file.path(FIG_DIR, "consensus_umap.pdf"), p_umap, width = 8, height = 6)
cat(sprintf("      Saved → figures/03_consensus/consensus_umap.pdf\n"))

cat("\n[done] Phase 3 complete.\n")
cat(sprintf("       High-confidence arm losses: %s\n",
            paste(arm_loss_freq$lost_arm, collapse = ", ")))
cat("\nNext: Rscript pipeline/04_expression_validation/expression_vs_cnv.R\n")
