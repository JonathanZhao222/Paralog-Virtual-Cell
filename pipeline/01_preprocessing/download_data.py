#!/usr/bin/env python3
"""
Download single-cell RNA-seq cancer data and reference files
for aneuploidy / paralog dependency analysis.

Sources
-------
- scRNA-seq      : CZ CellxGene Census  (cellxgene_census Python package)
- Gene positions : UCSC hg38 RefSeq     (required by InferCNV)
- Chromosome arms: UCSC hg38 cytoBand   (for arm-level CNV reporting)

Usage
-----
    python download_data.py                            # breast cancer, 30k cells
    python download_data.py --tumor "lung adenocarcinoma" --max_cells 20000
    python download_data.py --list_tumors              # show available options
    python download_data.py --skip_scrna               # reference files only
"""

import gzip
import os
import argparse
import requests
import pandas as pd
import scanpy as sc

# ---------------------------------------------------------------------------
# Paths  (all relative to project root, two levels up from this script)
# ---------------------------------------------------------------------------
ROOT      = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
SCRNASEQ  = os.path.join(ROOT, "data", "raw", "scrnaseq")
REFERENCE = os.path.join(ROOT, "data", "reference")
ENSEMBL   = os.path.join(ROOT, "data", "ensembl")
PROCESSED = os.path.join(ROOT, "data", "processed")

# ---------------------------------------------------------------------------
# Supported tumour types  (CellxGene Census disease labels)
# ---------------------------------------------------------------------------
SUPPORTED_TUMORS = [
    "breast cancer",
    "lung adenocarcinoma",
    "colorectal cancer",
    "ovarian cancer",
    "glioblastoma",
    "hepatocellular carcinoma",
    "bladder urothelial carcinoma",
    "pancreatic ductal adenocarcinoma",
]

# ---------------------------------------------------------------------------
# Reference file URLs  (hg38)
# ---------------------------------------------------------------------------
CYTOBAND_URL = "https://hgdownload.soe.ucsc.edu/goldenPath/hg38/database/cytoBand.txt.gz"
REFGENE_URL  = "https://hgdownload.soe.ucsc.edu/goldenPath/hg38/database/ncbiRefSeq.txt.gz"


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

def create_dirs():
    for d in [SCRNASEQ, REFERENCE, ENSEMBL, PROCESSED]:
        os.makedirs(d, exist_ok=True)
    print("[setup] Directory structure ready.")


def _download_gz(url: str, out_path: str, label: str):
    """Download a .gz file, decompress it, and save to out_path."""
    if os.path.exists(out_path):
        print(f"[{label}] Already exists — skipping ({out_path})")
        return

    print(f"[{label}] Downloading {url} ...")
    resp = requests.get(url, stream=True, timeout=180)
    resp.raise_for_status()

    gz_tmp = out_path + ".gz"
    with open(gz_tmp, "wb") as fh:
        for chunk in resp.iter_content(chunk_size=1 << 20):
            fh.write(chunk)

    print(f"[{label}] Decompressing ...")
    with gzip.open(gz_tmp, "rb") as gz_in, open(out_path, "wb") as f_out:
        f_out.write(gz_in.read())
    os.remove(gz_tmp)
    print(f"[{label}] Saved → {out_path}")


# ---------------------------------------------------------------------------
# 1. scRNA-seq  (CellxGene Census)
# ---------------------------------------------------------------------------

def download_scrna(tumor_type: str, max_cells: int):
    """
    Query CellxGene Census for primary cells from the given tumour type
    and save the result as a compressed AnnData (.h5ad) file.

    The raw count matrix is kept intact — normalisation happens in the
    Seurat QC script.
    """
    try:
        import cellxgene_census
    except ImportError:
        raise ImportError(
            "cellxgene-census is not installed.\n"
            "Run:  pip install cellxgene-census"
        )

    safe_name = tumor_type.replace(" ", "_")
    out_path  = os.path.join(SCRNASEQ, f"{safe_name}_raw.h5ad")

    if os.path.exists(out_path):
        print(f"[scRNA-seq] Already downloaded — skipping ({out_path})")
        return

    print(f"[scRNA-seq] Opening CellxGene Census ...")
    with cellxgene_census.open_soma() as census:

        # First: check how many cells are available
        obs_df = census["census_data"]["homo_sapiens"].obs.read(
            value_filter=f"disease == '{tumor_type}' and is_primary_data == True",
            column_names=["soma_joinid", "cell_type", "tissue",
                          "disease", "dataset_id", "donor_id", "assay"],
        ).concat().to_pandas()

        n_available = len(obs_df)
        if n_available == 0:
            print(f"[scRNA-seq] No cells found for '{tumor_type}'.")
            print("  Available tumour types:  " + ",  ".join(SUPPORTED_TUMORS))
            print("  Or run with --list_tumors to query Census directly.")
            return

        print(f"[scRNA-seq] Found {n_available:,} cells for '{tumor_type}'.")

        # Subsample if needed
        if n_available > max_cells:
            print(f"[scRNA-seq] Subsampling to {max_cells:,} cells ...")
            obs_df = obs_df.sample(n=max_cells, random_state=42)

        coords = obs_df["soma_joinid"].tolist()

        print(f"[scRNA-seq] Fetching count matrix for {len(coords):,} cells ...")
        adata = cellxgene_census.get_anndata(
            census,
            organism="Homo sapiens",
            obs_coords=coords,
        )

    print(f"[scRNA-seq] {adata.n_obs:,} cells × {adata.n_vars:,} genes")
    print(f"[scRNA-seq] Saving → {out_path}")
    adata.write_h5ad(out_path)
    print("[scRNA-seq] Done.\n")

    # Print a summary of datasets included
    datasets = obs_df["dataset_id"].value_counts()
    print("[scRNA-seq] Datasets included:")
    for ds_id, count in datasets.items():
        print(f"  {ds_id}  ({count:,} cells)")
    print()


# ---------------------------------------------------------------------------
# 2. Chromosome arm coordinates  (UCSC cytoBand → BED)
# ---------------------------------------------------------------------------

def download_chromosome_arms():
    """
    Download UCSC cytoBand (hg38) and convert to a simple BED file with
    one row per chromosome arm (e.g. chr1p, chr1q, ...).

    Output: data/reference/chromosome_arms_hg38.bed
    """
    raw_path = os.path.join(REFERENCE, "cytoBand_hg38.txt")
    arm_path = os.path.join(REFERENCE, "chromosome_arms_hg38.bed")

    _download_gz(CYTOBAND_URL, raw_path, "cytoBand")

    if os.path.exists(arm_path):
        print(f"[cytoBand] Arm BED already exists — skipping ({arm_path})")
        return

    df = pd.read_csv(
        raw_path, sep="\t", header=None,
        names=["chrom", "start", "end", "band", "stain"],
    )

    # Keep standard chromosomes only
    df = df[df["chrom"].str.match(r"^chr(\d+|X|Y)$")]
    df["arm"] = df["band"].str[0]   # 'p' or 'q'

    arms = (
        df.groupby(["chrom", "arm"])
          .agg(start=("start", "min"), end=("end", "max"))
          .reset_index()
    )
    arms["name"] = arms["chrom"] + arms["arm"]
    arms[["chrom", "start", "end", "name"]].to_csv(
        arm_path, sep="\t", index=False, header=False
    )
    print(f"[cytoBand] Chromosome arm BED saved → {arm_path}  ({len(arms)} arms)\n")


# ---------------------------------------------------------------------------
# 3. Gene positions  (UCSC RefSeq → InferCNV gene-order file)
# ---------------------------------------------------------------------------

def download_gene_positions():
    """
    Download UCSC ncbiRefSeq (hg38) and reformat to the 4-column
    tab-separated format required by InferCNV:

        gene_name   chrom   start   end

    One entry per gene (longest transcript kept).
    Output: data/reference/gene_order_hg38.txt
    """
    raw_path = os.path.join(REFERENCE, "ncbiRefSeq_hg38.txt")
    out_path = os.path.join(REFERENCE, "gene_order_hg38.txt")

    _download_gz(REFGENE_URL, raw_path, "genePos")

    if os.path.exists(out_path):
        print(f"[genePos] Gene order file already exists — skipping ({out_path})")
        return

    # UCSC ncbiRefSeq columns (no header in file):
    cols = [
        "bin", "name", "chrom", "strand",
        "txStart", "txEnd", "cdsStart", "cdsEnd",
        "exonCount", "exonStarts", "exonEnds",
        "score", "name2", "cdsStartStat", "cdsEndStat", "exonFrames",
    ]
    df = pd.read_csv(raw_path, sep="\t", header=None, names=cols, low_memory=False)

    # Keep standard chromosomes
    df = df[df["chrom"].str.match(r"^chr(\d+|X|Y)$")]

    df = df[["name2", "chrom", "txStart", "txEnd"]].rename(
        columns={"name2": "gene", "txStart": "start", "txEnd": "end"}
    )

    # One row per gene symbol — keep longest transcript
    df["length"] = df["end"] - df["start"]
    df = (
        df.sort_values("length", ascending=False)
          .drop_duplicates("gene")
          .sort_values(["chrom", "start"])
    )

    df[["gene", "chrom", "start", "end"]].to_csv(
        out_path, sep="\t", index=False, header=False
    )
    print(f"[genePos] InferCNV gene order file saved → {out_path}  ({len(df):,} genes)\n")


# ---------------------------------------------------------------------------
# 4. List available tumour types in Census  (optional helper)
# ---------------------------------------------------------------------------

def list_census_tumors():
    """Print all disease labels in Census that have > 1000 cells."""
    try:
        import cellxgene_census
    except ImportError:
        raise ImportError("Run:  pip install cellxgene-census")

    print("[Census] Querying available cancer types (this may take a minute) ...")
    with cellxgene_census.open_soma() as census:
        obs = census["census_data"]["homo_sapiens"].obs.read(
            column_names=["disease"]
        ).concat().to_pandas()

    counts = obs["disease"].value_counts()
    cancer_terms = counts[
        counts.index.str.contains("cancer|carcinoma|melanoma|glioma|lymphoma|leukemia",
                                  case=False, na=False)
    ]
    print(f"\n{'Disease':<50}  {'Cells':>10}")
    print("-" * 62)
    for disease, n in cancer_terms.items():
        print(f"{disease:<50}  {n:>10,}")


# ---------------------------------------------------------------------------
# Entry point
# ---------------------------------------------------------------------------

def main():
    parser = argparse.ArgumentParser(
        description="Download scRNA-seq + reference files for aneuploidy analysis",
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    parser.add_argument(
        "--tumor",
        default="breast cancer",
        help=(
            f"Tumour type to download from CellxGene Census. "
            f"Default: 'breast cancer'. "
            f"Use --list_tumors to see all options."
        ),
    )
    parser.add_argument(
        "--max_cells",
        type=int,
        default=30000,
        help=(
            "Max number of cells to download (default: 30,000). "
            "Larger = more accurate CNV calls but slower tools. "
            "Recommended range: 10,000–50,000."
        ),
    )
    parser.add_argument(
        "--skip_scrna",
        action="store_true",
        help="Skip scRNA-seq download and only download reference files.",
    )
    parser.add_argument(
        "--list_tumors",
        action="store_true",
        help="List all cancer types available in CellxGene Census and exit.",
    )
    args = parser.parse_args()

    if args.list_tumors:
        list_census_tumors()
        return

    create_dirs()

    if not args.skip_scrna:
        download_scrna(args.tumor, args.max_cells)

    download_chromosome_arms()
    download_gene_positions()

    print("=" * 60)
    print("Download complete. Files written to:")
    print(f"  data/raw/scrnaseq/   — count matrix (.h5ad)")
    print(f"  data/reference/      — gene order + chromosome arms")
    print()
    print("Next step:")
    print("  Rscript pipeline/01_preprocessing/seurat_qc_cluster.R")
    print("=" * 60)


if __name__ == "__main__":
    main()
