# Phase 1 — Seurat QC, normalisation, and clustering
#
# Input:   data/processed/10x/           (from convert_to_10x.py)
#          data/processed/cell_metadata.csv
# Output:  data/processed/seurat_qc.rds  (filtered + clustered Seurat object)
#          figures/01_preprocessing/     (QC and UMAP plots)

library(Seurat)
library(ggplot2)
library(dplyr)
library(patchwork)

# ── Paths ─────────────────────────────────────────────────────────────────────
ROOT     <- here::here()   # project root
TEN_X    <- file.path(ROOT, "data", "processed", "10x")
META_CSV <- file.path(ROOT, "data", "processed", "cell_metadata.csv")
OUT_RDS  <- file.path(ROOT, "data", "processed", "seurat_qc.rds")
FIG_DIR  <- file.path(ROOT, "figures", "01_preprocessing")
dir.create(FIG_DIR, showWarnings = FALSE, recursive = TRUE)

# ── 1. Load data ──────────────────────────────────────────────────────────────
cat("[1/7] Loading 10X data...\n")
counts <- Read10X(TEN_X)
seurat <- CreateSeuratObject(counts = counts, project = "breast_cancer",
                             min.cells = 3, min.features = 200)
cat(sprintf("      %d cells × %d genes loaded\n", ncol(seurat), nrow(seurat)))

# Add metadata from CellxGene (cell type, tissue, donor, etc.)
meta <- read.csv(META_CSV, row.names = 1)
shared_cells <- intersect(rownames(meta), colnames(seurat))
seurat <- AddMetaData(seurat, meta[shared_cells, ])

# ── 2. QC metrics ─────────────────────────────────────────────────────────────
cat("[2/7] Computing QC metrics...\n")
seurat[["percent.mt"]] <- PercentageFeatureSet(seurat, pattern = "^MT-")

p_qc <- VlnPlot(seurat,
                features = c("nFeature_RNA", "nCount_RNA", "percent.mt"),
                ncol = 3, pt.size = 0) &
        theme(axis.text.x = element_blank())
ggsave(file.path(FIG_DIR, "qc_violin_prefilter.pdf"), p_qc,
       width = 12, height = 4)

p_scatter <- FeatureScatter(seurat, "nCount_RNA", "nFeature_RNA") +
             geom_hline(yintercept = c(200, 6000), linetype = "dashed", colour = "red")
ggsave(file.path(FIG_DIR, "qc_count_vs_features.pdf"), p_scatter,
       width = 6, height = 5)

# ── 3. Filter ─────────────────────────────────────────────────────────────────
cat("[3/7] Filtering cells...\n")
n_before <- ncol(seurat)
seurat <- subset(seurat,
                 nFeature_RNA > 200 &
                 nFeature_RNA < 6000 &
                 percent.mt   < 25)
cat(sprintf("      %d → %d cells after QC filter\n", n_before, ncol(seurat)))

# ── 4. Normalise + variable features ─────────────────────────────────────────
cat("[4/7] Normalising and finding variable features...\n")
seurat <- NormalizeData(seurat, normalization.method = "LogNormalize",
                        scale.factor = 10000, verbose = FALSE)
seurat <- FindVariableFeatures(seurat, nfeatures = 2000, verbose = FALSE)

p_vf <- VariableFeaturePlot(seurat) +
        ggtitle("Top 2000 Variable Features")
ggsave(file.path(FIG_DIR, "variable_features.pdf"), p_vf, width = 8, height = 5)

# ── 5. Scale + PCA ────────────────────────────────────────────────────────────
cat("[5/7] Scaling and running PCA...\n")
seurat <- ScaleData(seurat, verbose = FALSE)
seurat <- RunPCA(seurat, npcs = 30, verbose = FALSE)

p_elbow <- ElbowPlot(seurat, ndims = 30) + ggtitle("PCA Elbow Plot")
ggsave(file.path(FIG_DIR, "pca_elbow.pdf"), p_elbow, width = 6, height = 4)

# ── 6. UMAP + clustering ──────────────────────────────────────────────────────
cat("[6/7] Running UMAP and clustering...\n")
seurat <- RunUMAP(seurat,   dims = 1:20, verbose = FALSE)
seurat <- FindNeighbors(seurat, dims = 1:20, verbose = FALSE)
seurat <- FindClusters(seurat,  resolution = 0.5, verbose = FALSE)

# Figure A: clusters
p_cluster <- DimPlot(seurat, reduction = "umap", group.by = "seurat_clusters",
                     label = TRUE, repel = TRUE) +
             ggtitle("Seurat Clusters") +
             theme_classic()

# Figure B: cell type from CellxGene metadata
p_celltype <- if ("cell_type" %in% colnames(seurat@meta.data)) {
  DimPlot(seurat, reduction = "umap", group.by = "cell_type",
          label = TRUE, repel = TRUE, label.size = 3) +
  ggtitle("Cell Type (CellxGene annotation)") +
  theme_classic() +
  theme(legend.position = "none")
} else {
  ggplot() + annotate("text", x = 0.5, y = 0.5, label = "No cell_type metadata") +
  theme_void()
}

# Figure C: mitochondrial % on UMAP (quality check)
p_mt <- FeaturePlot(seurat, features = "percent.mt", reduction = "umap") +
        ggtitle("% Mitochondrial Reads") +
        scale_colour_viridis_c()

p_combined <- (p_cluster | p_celltype) / p_mt
ggsave(file.path(FIG_DIR, "umap_overview.pdf"), p_combined,
       width = 14, height = 12)

# ── 7. Label tumour vs. normal ────────────────────────────────────────────────
# CellxGene annotations separate epithelial/tumour from immune/stromal cells.
# Tumour cells are needed as the "test" group for CNV tools;
# immune/stromal cells serve as the normal reference.
cat("[7/7] Annotating tumour vs. normal cells...\n")

tumour_keywords <- c("epithelial", "luminal", "basal", "cancer",
                     "malignant", "tumor", "tumour")
if ("cell_type" %in% colnames(seurat@meta.data)) {
  seurat$tumour_status <- ifelse(
    grepl(paste(tumour_keywords, collapse = "|"),
          seurat$cell_type, ignore.case = TRUE),
    "tumour", "normal_reference"
  )
} else {
  # Fallback: flag as unknown — adjust manually after inspecting the UMAP
  seurat$tumour_status <- "unknown"
  warning("No cell_type metadata found. Set seurat$tumour_status manually.")
}

p_tumour <- DimPlot(seurat, reduction = "umap", group.by = "tumour_status",
                    cols = c("tumour" = "#E74C3C", "normal_reference" = "#3498DB",
                             "unknown" = "grey80")) +
            ggtitle("Tumour vs. Normal Reference") +
            theme_classic()
ggsave(file.path(FIG_DIR, "umap_tumour_vs_normal.pdf"), p_tumour,
       width = 7, height = 6)

tumour_n  <- sum(seurat$tumour_status == "tumour")
normal_n  <- sum(seurat$tumour_status == "normal_reference")
cat(sprintf("      Tumour cells: %d | Normal reference cells: %d\n",
            tumour_n, normal_n))

# ── Save ──────────────────────────────────────────────────────────────────────
saveRDS(seurat, OUT_RDS)
cat(sprintf("\n[done] Seurat object saved → %s\n", OUT_RDS))
cat("[done] Figures saved → figures/01_preprocessing/\n")
cat("\nNext step:\n  Rscript pipeline/02_cnv/run_infercnv.R\n")
