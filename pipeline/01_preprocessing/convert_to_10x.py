#!/usr/bin/env python3
"""
Convert the downloaded .h5ad file to 10X sparse format
so Seurat can read it natively with Read10X().

Also saves cell metadata as a separate CSV for adding back
into the Seurat object (cell type, tissue, donor, dataset).

Output
------
data/processed/10x/
    matrix.mtx.gz
    barcodes.tsv.gz
    features.tsv.gz
data/processed/cell_metadata.csv
"""

import os
import gzip
import anndata
import scipy.io
import pandas as pd
import numpy as np

ROOT      = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
H5AD_PATH = os.path.join(ROOT, "data", "raw", "scrnaseq", "breast_cancer_raw.h5ad")
OUT_DIR   = os.path.join(ROOT, "data", "processed", "10x")
META_PATH = os.path.join(ROOT, "data", "processed", "cell_metadata.csv")


def convert():
    os.makedirs(OUT_DIR, exist_ok=True)

    print(f"[convert] Loading {H5AD_PATH} ...")
    adata = anndata.read_h5ad(H5AD_PATH)
    print(f"[convert] {adata.n_obs:,} cells × {adata.n_vars:,} genes")

    # ── Filter to primary breast tissue only ──────────────────────────────────
    # The Census "breast cancer" disease label includes metastatic sites
    # (liver, brain, bone, lung). Keep breast + axilla (regional lymph nodes)
    # which represent the primary tumour and its immediate drainage.
    if "tissue" in adata.obs.columns:
        keep_tissues = ["breast", "axilla"]
        mask = adata.obs["tissue"].isin(keep_tissues)
        n_before = adata.n_obs
        adata = adata[mask].copy()
        print(f"[convert] Tissue filter: {n_before:,} → {adata.n_obs:,} cells "
              f"(kept: {', '.join(keep_tissues)})")
        print(f"[convert] Tissue breakdown after filter:")
        for tissue, n in adata.obs["tissue"].value_counts().items():
            print(f"  {tissue:<30} {n:>6,} cells")

    # ── 1. Count matrix ───────────────────────────────────────────────────────
    # Seurat expects genes × cells (transposed relative to AnnData)
    X = adata.X
    if not hasattr(X, "toarray"):
        from scipy.sparse import csr_matrix
        X = csr_matrix(X)
    X_T = X.T.tocsr()

    mtx_path = os.path.join(OUT_DIR, "matrix.mtx.gz")
    print(f"[convert] Writing matrix ({X_T.shape[0]} genes × {X_T.shape[1]} cells) ...")
    with gzip.open(mtx_path, "wb") as f:
        scipy.io.mmwrite(f, X_T)
    print(f"[convert] Saved → {mtx_path}")

    # ── 2. Barcodes (cell IDs) ────────────────────────────────────────────────
    bc_path = os.path.join(OUT_DIR, "barcodes.tsv.gz")
    with gzip.open(bc_path, "wt") as f:
        f.write("\n".join(adata.obs_names) + "\n")
    print(f"[convert] Saved → {bc_path}")

    # ── 3. Features (genes) ───────────────────────────────────────────────────
    # 10X format: gene_id \t gene_name \t feature_type
    feat_path = os.path.join(OUT_DIR, "features.tsv.gz")
    gene_ids   = adata.var.get("gene_ids",   adata.var_names)
    gene_names = adata.var.get("feature_name", adata.var_names)
    with gzip.open(feat_path, "wt") as f:
        for gid, gname in zip(gene_ids, gene_names):
            f.write(f"{gid}\t{gname}\tGene Expression\n")
    print(f"[convert] Saved → {feat_path}")

    # ── 4. Cell metadata CSV ──────────────────────────────────────────────────
    meta_cols = [c for c in ["cell_type", "tissue", "disease", "donor_id",
                              "dataset_id", "assay"] if c in adata.obs.columns]
    meta = adata.obs[meta_cols].copy()

    # Clean up cell type labels for use in R (no special characters)
    if "cell_type" in meta.columns:
        meta["cell_type_clean"] = (
            meta["cell_type"]
            .str.replace(r"[^A-Za-z0-9_]", "_", regex=True)
            .str.strip("_")
        )

    meta.to_csv(META_PATH)
    print(f"[convert] Metadata saved → {META_PATH}")
    print(f"[convert] Cell types present:")
    if "cell_type" in meta.columns:
        for ct, n in meta["cell_type"].value_counts().items():
            print(f"  {ct:<45} {n:>6,} cells")

    print("\n[convert] Done. Next step:")
    print("  Rscript pipeline/01_preprocessing/seurat_qc_cluster.R")


if __name__ == "__main__":
    convert()
