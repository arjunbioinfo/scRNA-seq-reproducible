#!/usr/bin/env python3
"""
Reproducible single-cell RNA-seq pipeline (Scanpy / scverse).

One config file drives the whole run: point it at a GEO accession, a local
10x/h5/h5ad file, a CSV count matrix, or a URL, and it performs
QC -> normalization -> HVG -> PCA ->
neighbors -> UMAP -> Leiden clustering -> marker genes, saving a processed
AnnData, figures, and CSV summaries.

Usage:
    python scrna_pipeline.py --config config/config.yaml

The pipeline is deterministic given the config (fixed random seeds).
"""
from __future__ import annotations
import argparse
import gzip
import os
import sys
import shutil
import tarfile
import urllib.request
from pathlib import Path

import numpy as np
import pandas as pd
import scanpy as sc
import yaml


# --------------------------------------------------------------------------- #
# Config
# --------------------------------------------------------------------------- #
def load_config(path: str) -> dict:
    with open(path) as fh:
        return yaml.safe_load(fh)


# --------------------------------------------------------------------------- #
# Data loading (GEO accession | local 10x/h5/h5ad | URL)
# --------------------------------------------------------------------------- #
def _download(url: str, dest: Path) -> Path:
    dest.parent.mkdir(parents=True, exist_ok=True)
    if not dest.exists():
        print(f"[io] downloading {url}")
        urllib.request.urlretrieve(url, dest)
    return dest


def _read_any(path: Path, var_names: str) -> sc.AnnData:
    """Read a single file by extension: .h5ad, 10x .h5, or a 10x mtx directory."""
    p = str(path)
    if p.endswith(".h5ad"):
        return sc.read_h5ad(p)
    if p.endswith(".h5"):
        return sc.read_10x_h5(p)
    if path.is_dir():
        return sc.read_10x_mtx(p, var_names=var_names, cache=False)
    raise ValueError(f"Unrecognized data path: {path}")


def load_data(cfg: dict, work_dir: Path) -> sc.AnnData:
    d = cfg["data"]
    source = d["source"]
    var_names = d.get("var_names", "gene_symbols")

    if source == "example":
        print("[io] loading example dataset: 3k PBMCs (10x)")
        adata = sc.datasets.pbmc3k()

    elif source == "h5ad":
        adata = sc.read_h5ad(d["path"])

    elif source in ("10x_mtx", "h5"):
        adata = _read_any(Path(d["path"]), var_names)

    elif source == "url":
        url = d["url"]
        fname = url.split("?")[0].split("/")[-1] or "download.h5ad"
        local = _download(url, work_dir / fname)
        adata = _read_any(local, var_names)

    elif source == "csv":
        adata = load_csv(d, work_dir)

    elif source == "geo":
        adata = load_geo(d["geo_accession"], work_dir, var_names)

    else:
        raise ValueError(f"Unknown data.source: {source!r}")

    adata.var_names_make_unique()
    adata.obs_names_make_unique()
    print(f"[io] loaded AnnData: {adata.n_obs} cells x {adata.n_vars} genes")
    return adata


def load_csv(d: dict, work_dir: Path) -> sc.AnnData:
    """Load a dense CSV/TSV count matrix (optionally .gz) into AnnData.

    Many GEO series ship counts as a single delimited table. By convention rows
    are genes and columns are cells (``genes_are_rows: true``); flip it to false
    if your file is cells x genes. ``path`` reads a local file, ``url`` downloads
    one first. Compression is inferred from a ``.gz`` suffix.
    """
    import scipy.sparse as sp

    src = d.get("path")
    if d.get("url"):
        url = d["url"]
        fname = url.split("?")[0].split("/")[-1] or "matrix.csv.gz"
        src = str(_download(url, work_dir / fname))
    if not src:
        raise ValueError("data.source 'csv' needs a 'path' or 'url'")

    sep = d.get("sep", ",")
    print(f"[io] reading CSV count matrix {src}")
    df = pd.read_csv(src, index_col=0, sep=sep)   # .gz auto-detected by pandas
    if d.get("genes_are_rows", True):
        df = df.T                                 # -> cells (rows) x genes (cols)
    adata = sc.AnnData(
        X=sp.csr_matrix(df.to_numpy(dtype="float32")),
        obs=pd.DataFrame(index=df.index.astype(str)),
        var=pd.DataFrame(index=df.columns.astype(str)),
    )
    print(f"[io] CSV -> {adata.n_obs} cells x {adata.n_vars} genes")
    return adata


def load_geo(accession: str, work_dir: Path, var_names: str) -> sc.AnnData:
    """
    Download a GEO series' supplementary files and auto-detect a loadable matrix.
    GEO supplementary layouts vary; this handles the common cases:
      - a .h5ad or 10x .h5 supplementary file
      - a 10x mtx triplet (matrix.mtx(.gz), (features|genes).tsv(.gz), barcodes.tsv(.gz))
      - a .tar of the above
    If auto-detection fails, download manually and use source: 10x_mtx / h5 / h5ad.
    """
    try:
        import GEOparse  # noqa
    except Exception as exc:  # pragma: no cover
        raise RuntimeError(
            "GEO loading needs GEOparse (in environment.yml). "
            f"Import failed: {exc}"
        )
    import GEOparse

    dest = work_dir / accession
    dest.mkdir(parents=True, exist_ok=True)
    print(f"[io] fetching GEO {accession} (supplementary files)")
    gse = GEOparse.get_GEO(geo=accession, destdir=str(dest), how="brief", silent=True)
    gse.download_supplementary_files(directory=str(dest), download_sra=False)

    files = list(dest.rglob("*"))
    # 1) direct h5ad / 10x h5
    for f in files:
        if f.suffix == ".h5ad" or f.name.endswith(".h5"):
            return _read_any(f, var_names)
    # 2) untar any tarballs
    for f in files:
        if f.name.endswith((".tar", ".tar.gz", ".tgz")):
            with tarfile.open(f) as t:
                t.extractall(dest)
    # 3) locate a 10x mtx directory (a folder containing matrix.mtx[.gz])
    for mtx in dest.rglob("matrix.mtx*"):
        folder = mtx.parent
        # scanpy expects features/genes + barcodes alongside; normalize names
        return sc.read_10x_mtx(str(folder), var_names=var_names, cache=False)

    raise RuntimeError(
        f"Could not auto-detect a matrix in {accession}. Download the supplementary "
        f"file manually and set data.source to 10x_mtx / h5 / h5ad in the config."
    )


# --------------------------------------------------------------------------- #
# QC  (scverse best practices: MAD-based outlier filtering)
# --------------------------------------------------------------------------- #
def _is_outlier(adata, metric: str, nmads: int) -> pd.Series:
    M = adata.obs[metric]
    med = np.median(M)
    mad = np.median(np.abs(M - med))
    return (M < med - nmads * mad) | (M > med + nmads * mad)


def run_qc(adata, cfg: dict, figdir: Path):
    q = cfg["qc"]
    adata.var["mt"] = adata.var_names.str.startswith(q.get("mito_prefix", "MT-"))
    adata.var["ribo"] = adata.var_names.str.startswith(("RPS", "RPL"))
    adata.var["hb"] = adata.var_names.str.contains(r"^HB[^(P)]", regex=True)
    sc.pp.calculate_qc_metrics(
        adata, qc_vars=["mt", "ribo", "hb"], inplace=True, percent_top=[20], log1p=True
    )

    sc.pl.violin(
        adata, ["n_genes_by_counts", "total_counts", "pct_counts_mt"],
        jitter=0.4, multi_panel=True, show=False, save="_qc.png",
    )
    _move_fig(figdir, "violin_qc.png")

    n_before = adata.n_obs
    nmads = q.get("n_mads", 5)
    adata.obs["outlier"] = (
        _is_outlier(adata, "log1p_total_counts", nmads)
        | _is_outlier(adata, "log1p_n_genes_by_counts", nmads)
        | _is_outlier(adata, "pct_counts_in_top_20_genes", nmads)
    )
    mt_cap = q.get("pct_counts_mt_max")
    if mt_cap is not None:
        adata.obs["mt_outlier"] = adata.obs["pct_counts_mt"] > mt_cap
    else:
        adata.obs["mt_outlier"] = _is_outlier(adata, "pct_counts_mt", 3)

    adata = adata[~(adata.obs["outlier"] | adata.obs["mt_outlier"])].copy()
    sc.pp.filter_cells(adata, min_genes=q.get("min_genes", 200))
    sc.pp.filter_genes(adata, min_cells=q.get("min_cells", 3))
    print(f"[qc] cells: {n_before} -> {adata.n_obs} (removed {n_before - adata.n_obs})")

    if q.get("doublet_detection", True):
        try:
            sc.pp.scrublet(adata, random_state=0)
            adata = adata[~adata.obs["predicted_doublet"]].copy()
            print(f"[qc] after doublet removal: {adata.n_obs} cells")
        except Exception as exc:
            print(f"[qc] scrublet skipped ({exc})")
    return adata


# --------------------------------------------------------------------------- #
# Normalization / HVG / PCA
# --------------------------------------------------------------------------- #
def preprocess(adata, cfg: dict, figdir: Path):
    p = cfg["preprocess"]
    adata.layers["counts"] = adata.X.copy()
    sc.pp.normalize_total(adata, target_sum=p.get("target_sum", 1e4))
    sc.pp.log1p(adata)
    adata.raw = adata

    sc.pp.highly_variable_genes(adata, n_top_genes=p.get("n_top_genes", 2000))
    sc.pl.highly_variable_genes(adata, show=False, save="_hvg.png")
    _move_fig(figdir, "hvg.png")
    adata = adata[:, adata.var["highly_variable"]].copy()

    for cov in p.get("regress_out", []) or []:
        pass
    if p.get("regress_out"):
        sc.pp.regress_out(adata, p["regress_out"])
    sc.pp.scale(adata, max_value=p.get("scale_max", 10))
    sc.tl.pca(adata, n_comps=p.get("n_pcs", 50), svd_solver="arpack", random_state=0)
    sc.pl.pca_variance_ratio(adata, n_pcs=p.get("n_pcs", 50), show=False, save="_pca.png")
    _move_fig(figdir, "pca_variance.png")
    return adata


# --------------------------------------------------------------------------- #
# Neighbors / UMAP / Leiden / markers
# --------------------------------------------------------------------------- #
def cluster(adata, cfg: dict, figdir: Path):
    c = cfg["cluster"]
    sc.pp.neighbors(
        adata, n_neighbors=c.get("n_neighbors", 15),
        n_pcs=cfg["preprocess"].get("n_pcs", 50), random_state=c.get("random_state", 0),
    )
    sc.tl.umap(adata, random_state=c.get("random_state", 0))
    sc.tl.leiden(
        adata, resolution=c.get("leiden_resolution", 1.0),
        random_state=c.get("random_state", 0), flavor="igraph", n_iterations=2, directed=False,
    )
    print(f"[cluster] {adata.obs['leiden'].nunique()} Leiden clusters")
    sc.pl.umap(adata, color=["leiden"], legend_loc="on data", show=False, save="_leiden.png")
    _move_fig(figdir, "umap_leiden.png")
    sc.pl.umap(
        adata, color=["total_counts", "pct_counts_mt"], show=False, save="_qc_umap.png"
    )
    _move_fig(figdir, "umap_qc.png")
    return adata


def find_markers(adata, cfg: dict, figdir: Path, outdir: Path):
    m = cfg["markers"]
    sc.tl.rank_genes_groups(
        adata, "leiden", method=m.get("method", "wilcoxon"), use_raw=True
    )
    sc.pl.rank_genes_groups_dotplot(
        adata, n_genes=min(5, m.get("n_genes", 25)), show=False, save="_markers.png"
    )
    _move_fig(figdir, "markers_dotplot.png")
    df = sc.get.rank_genes_groups_df(adata, group=None)
    df.to_csv(outdir / "markers.csv", index=False)
    (adata.obs["leiden"].value_counts().sort_index()
     .rename_axis("cluster").rename("n_cells")
     .to_csv(outdir / "cluster_sizes.csv"))
    print(f"[markers] wrote {outdir/'markers.csv'} and cluster_sizes.csv")
    return adata


# --------------------------------------------------------------------------- #
# helpers
# --------------------------------------------------------------------------- #
def _move_fig(figdir: Path, name: str):
    """Scanpy writes into ./figures/; move the newest file to our figures dir."""
    src_dir = Path("figures")
    if not src_dir.exists():
        return
    figs = sorted(src_dir.glob("*"), key=os.path.getmtime)
    if figs:
        figdir.mkdir(parents=True, exist_ok=True)
        shutil.move(str(figs[-1]), str(figdir / name))


# --------------------------------------------------------------------------- #
# main
# --------------------------------------------------------------------------- #
def main():
    ap = argparse.ArgumentParser(description="Reproducible scRNA-seq pipeline (Scanpy)")
    ap.add_argument("--config", default="config/config.yaml")
    args = ap.parse_args()

    cfg = load_config(args.config)
    sc.settings.verbosity = 1
    np.random.seed(0)

    out = Path(cfg["output"]["dir"]); out.mkdir(parents=True, exist_ok=True)
    figdir = Path(cfg["output"]["figures_dir"]); figdir.mkdir(parents=True, exist_ok=True)
    work = Path(cfg["output"].get("work_dir", "data")); work.mkdir(parents=True, exist_ok=True)
    sc.settings.figdir = "figures"

    adata = load_data(cfg, work)
    adata = run_qc(adata, cfg, figdir)
    adata = preprocess(adata, cfg, figdir)
    adata = cluster(adata, cfg, figdir)
    adata = find_markers(adata, cfg, figdir, out)

    processed = Path(cfg["output"]["processed_h5ad"])
    processed.parent.mkdir(parents=True, exist_ok=True)
    adata.write_h5ad(processed)
    print(f"[done] processed AnnData -> {processed}")
    print(f"[done] figures -> {figdir}/  | tables -> {out}/")


if __name__ == "__main__":
    sys.exit(main())
