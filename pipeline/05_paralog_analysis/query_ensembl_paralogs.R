# Phase 5a — Paralog Query via Ensembl BioMart
#
# For each gene on a high-confidence lost arm, retrieve its human paralogs
# from Ensembl BioMart. Keep only cross-arm paralog pairs — i.e., the paralog
# resides on a DIFFERENT chromosome arm than the lost gene.
#
# The biological hypothesis: when arm X is lost in tumour cells, a gene G on X
# may be compensated by its paralog P on a different arm (which is intact). This
# creates a synthetic dependency: cells that have lost arm X become dependent on P.
#
# Input:
#   results/consensus/arm_loss_freq.csv     — arms with high-confidence losses
#   data/reference/gene_order_hg38.txt      — gene → chromosome position
#   data/reference/chromosome_arms_hg38.bed — arm boundary definitions
#
# Output:
#   results/paralog_pairs/all_paralogs.csv          — raw BioMart query results
#   results/paralog_pairs/cross_arm_paralogs.csv    — filtered: paralog on different arm
#   results/paralog_pairs/paralog_summary.csv       — per-lost-arm summary counts

library(biomaRt)
library(dplyr)

# ── Paths ──────────────────────────────────────────────────────────────────────
ROOT       <- getwd()
ARM_FREQ   <- file.path(ROOT, "results", "consensus", "arm_loss_freq.csv")
GENE_ORDER <- file.path(ROOT, "data", "reference", "gene_order_hg38.txt")
ARM_BED    <- file.path(ROOT, "data", "reference", "chromosome_arms_hg38.bed")
OUT_DIR    <- file.path(ROOT, "results", "paralog_pairs")
dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)

cat("=== Phase 5a: Ensembl Paralog Query ===\n\n")

# ── 1. Load arm-loss data ──────────────────────────────────────────────────────
cat("[1/4] Loading lost-arm gene lists...\n")
arm_freq   <- read.csv(ARM_FREQ)
gene_order <- read.table(GENE_ORDER, col.names = c("gene", "chrom", "start", "end"))
arm_bed    <- read.table(ARM_BED,    col.names = c("chrom", "start", "end", "arm"))

# Assign arm to each gene
assign_arm <- function(chrom, pos) {
  hits <- arm_bed[arm_bed$chrom == chrom & arm_bed$start <= pos & arm_bed$end > pos, ]
  if (nrow(hits) == 0) return(NA_character_)
  hits$arm[1]
}
gene_order$arm <- mapply(assign_arm, gene_order$chrom,
                         floor((gene_order$start + gene_order$end) / 2))
gene_order <- gene_order[!is.na(gene_order$arm), ]

# Genes on lost arms
lost_genes <- gene_order[gene_order$arm %in% arm_freq$lost_arm, ]
cat(sprintf("      %d genes on %d high-confidence lost arms\n",
            nrow(lost_genes), length(unique(lost_genes$arm))))

# ── 2. Connect to Ensembl BioMart ─────────────────────────────────────────────
cat("[2/4] Connecting to Ensembl BioMart (GRCh38)...\n")
mart <- tryCatch(
  useEnsembl("ensembl", dataset = "hsapiens_gene_ensembl"),
  error = function(e) {
    cat("      Primary mirror failed — trying US mirror...\n")
    tryCatch(
      useEnsembl("ensembl", dataset = "hsapiens_gene_ensembl", mirror = "useast"),
      error = function(e2) {
        cat("      US mirror failed — trying Asia mirror...\n")
        useEnsembl("ensembl", dataset = "hsapiens_gene_ensembl", mirror = "asia")
      }
    )
  }
)
cat("      Connected.\n")

# ── 3. Query paralogs ─────────────────────────────────────────────────────────
# BioMart does not allow mixing the "Gene" attributes page (hgnc_symbol) with
# the "Homologs" page (hsapiens_paralog_*) in a single query.
# Solution: two-step query.
#   Step A: hgnc_symbol → ensembl_gene_id  (Gene page)
#   Step B: ensembl_gene_id → paralog attrs (Homologs page)
#   Then join the two tables on ensembl_gene_id.

cat("[3/4] Querying paralogs (two-step BioMart approach)...\n")
cat("      This may take 2-10 minutes depending on Ensembl server load.\n")

gene_list  <- lost_genes$gene
BATCH_SIZE <- 500   # larger batches are fine for single-page queries

# ── Step A: HGNC symbol → Ensembl gene ID ────────────────────────────────────
cat("      Step A: mapping HGNC symbols to Ensembl IDs...\n")
n_batches_a <- ceiling(length(gene_list) / BATCH_SIZE)
id_maps <- list()
for (i in seq_len(n_batches_a)) {
  batch <- gene_list[((i-1)*BATCH_SIZE + 1):min(i*BATCH_SIZE, length(gene_list))]
  res <- tryCatch(
    getBM(attributes = c("hgnc_symbol", "ensembl_gene_id"),
          filters    = "hgnc_symbol",
          values     = batch,
          mart       = mart),
    error = function(e) { cat(sprintf("      WARNING A batch %d: %s\n", i, conditionMessage(e))); NULL }
  )
  if (!is.null(res) && nrow(res) > 0) id_maps[[i]] <- res
  Sys.sleep(0.3)
}
gene_map <- unique(do.call(rbind, Filter(Negate(is.null), id_maps)))
gene_map <- gene_map[gene_map$hgnc_symbol != "" & gene_map$ensembl_gene_id != "", ]
cat(sprintf("      Mapped %d / %d genes to Ensembl IDs\n",
            nrow(gene_map), length(gene_list)))

# ── Step B: Ensembl IDs → paralog attributes (Homologs page only) ────────────
cat("      Step B: querying paralog attributes...\n")
paralog_attrs <- c(
  "ensembl_gene_id",
  "hsapiens_paralog_associated_gene_name",
  "hsapiens_paralog_chromosome",
  "hsapiens_paralog_chrom_start",
  "hsapiens_paralog_chrom_end",
  "hsapiens_paralog_perc_id",
  "hsapiens_paralog_perc_id_r1",
  "hsapiens_paralog_orthology_type"
)

ensembl_ids <- gene_map$ensembl_gene_id
n_batches_b <- ceiling(length(ensembl_ids) / BATCH_SIZE)
all_paralogs <- list()
for (i in seq_len(n_batches_b)) {
  batch <- ensembl_ids[((i-1)*BATCH_SIZE + 1):min(i*BATCH_SIZE, length(ensembl_ids))]
  cat(sprintf("      Batch %d / %d (%d genes)...\n", i, n_batches_b, length(batch)))
  res <- tryCatch(
    getBM(attributes = paralog_attrs,
          filters    = "ensembl_gene_id",
          values     = batch,
          mart       = mart),
    error = function(e) { cat(sprintf("      WARNING B batch %d: %s\n", i, conditionMessage(e))); NULL }
  )
  if (!is.null(res) && nrow(res) > 0) all_paralogs[[i]] <- res
  Sys.sleep(0.3)
}

paralogs_raw <- do.call(rbind, Filter(Negate(is.null), all_paralogs))
colnames(paralogs_raw) <- c(
  "ensembl_gene_id", "paralog_gene", "paralog_chrom",
  "paralog_start", "paralog_end",
  "perc_id_fwd", "perc_id_rev",
  "paralog_type"
)

# Join back to get HGNC gene names
paralogs_raw <- merge(paralogs_raw, gene_map, by = "ensembl_gene_id", all.x = TRUE)
paralogs_raw <- paralogs_raw[, c("hgnc_symbol", "paralog_gene", "paralog_chrom",
                                  "paralog_start", "paralog_end",
                                  "perc_id_fwd", "perc_id_rev", "paralog_type")]
colnames(paralogs_raw)[1] <- "gene"

# Filter: remove self-hits and missing data
paralogs_raw <- paralogs_raw[
  !is.na(paralogs_raw$paralog_gene) &
  paralogs_raw$paralog_gene != "" &
  paralogs_raw$paralog_gene != paralogs_raw$gene,
]

cat(sprintf("      Raw paralog pairs: %d\n", nrow(paralogs_raw)))
write.csv(paralogs_raw, file.path(OUT_DIR, "all_paralogs.csv"), row.names = FALSE)
cat("      Saved → results/paralog_pairs/all_paralogs.csv\n")

# ── 4. Filter to cross-arm paralog pairs ──────────────────────────────────────
cat("[4/4] Filtering to cross-arm paralog pairs...\n")

# Assign arm to each query gene
paralogs_raw <- paralogs_raw %>%
  left_join(gene_order %>% select(gene, arm) %>% rename(query_arm = arm),
            by = "gene")

# Assign arm to each paralog gene
paralog_arm_map <- gene_order %>%
  select(gene, arm) %>%
  rename(paralog_gene = gene, paralog_arm = arm)
paralogs_raw <- paralogs_raw %>%
  left_join(paralog_arm_map, by = "paralog_gene")

# Cross-arm = query and paralog on different arms
cross_arm <- paralogs_raw %>%
  filter(
    !is.na(query_arm) & !is.na(paralog_arm) &
    query_arm != paralog_arm &
    query_arm %in% arm_freq$lost_arm   # query gene is on a lost arm
  ) %>%
  # Paralog arm must NOT also be a commonly-lost arm
  # (we want the paralog to be on an intact arm)
  filter(!paralog_arm %in% arm_freq$lost_arm)

cat(sprintf("      Cross-arm pairs (paralog on intact arm): %d\n", nrow(cross_arm)))
write.csv(cross_arm, file.path(OUT_DIR, "cross_arm_paralogs.csv"), row.names = FALSE)
cat("      Saved → results/paralog_pairs/cross_arm_paralogs.csv\n")

# ── Summary per lost arm ───────────────────────────────────────────────────────
summary_tbl <- cross_arm %>%
  group_by(query_arm) %>%
  summarise(
    n_genes_with_paralogs = n_distinct(gene),
    n_paralog_pairs       = n(),
    top_paralogs          = paste(head(unique(paralog_gene), 5), collapse = ", "),
    .groups = "drop"
  )
print(summary_tbl)
write.csv(summary_tbl, file.path(OUT_DIR, "paralog_summary.csv"), row.names = FALSE)
cat("\nSaved → results/paralog_pairs/paralog_summary.csv\n")
cat("\n[done] Phase 5a complete.\n")
cat("\nNext: Rscript pipeline/05_paralog_analysis/paralog_expression.R\n")
