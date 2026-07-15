# scRNA-seq-reproducible

[![CI](https://github.com/arjunbioinfo/scRNA-seq-reproducible/actions/workflows/ci.yml/badge.svg)](https://github.com/arjunbioinfo/scRNA-seq-reproducible/actions/workflows/ci.yml)

A reproducible, config-driven single-cell RNA-seq pipeline built on
[Scanpy](https://scanpy.readthedocs.io) / [scverse](https://scverse.org). Point it
at **any** public dataset — a GEO accession, a local 10x/`.h5`/`.h5ad` file, or a
URL — and it runs the standard workflow end to end and writes a processed
`AnnData`, figures, and marker tables.

Pipeline: **load → QC (MAD-based filtering + doublet removal) → normalize → HVG →
PCA → neighbors → UMAP → Leiden clustering → marker genes**. It is deterministic
given the config (fixed random seeds), and versions are pinned in
`environment.yml` / `requirements.txt`.

## Quick start

```bash
# 1. environment (conda)
conda env create -f environment.yml
conda activate scrna-repro
#    ...or pip:  pip install -r requirements.txt

# 2. run on the bundled example (3k PBMCs, downloads automatically)
python scrna_pipeline.py --config config/config.yaml
#    ...or:  make run
```

Outputs land in `results/` (processed `.h5ad`, `markers.csv`, `cluster_sizes.csv`)
and `figures/` (QC violins, HVG, PCA variance, UMAPs, marker dotplot).

Prefer an interactive walkthrough? Open **`notebooks/walkthrough.ipynb`** — the
same steps, cell by cell, with plots inline.

Every push runs the pipeline on the example dataset in CI
(`.github/workflows/ci.yml`), so the badge above shows whether the workflow is
reproducible from a clean environment.

## Use your own data

Edit the `data` block in `config/config.yaml`:

```yaml
# a) GEO accession (auto-download + auto-detect matrix)
data: { source: geo, geo_accession: GSE123456 }

# b) local 10x cellranger output
data: { source: 10x_mtx, path: data/filtered_feature_bc_matrix, var_names: gene_symbols }

# c) local 10x .h5 or an .h5ad
data: { source: h5,   path: data/sample.h5 }
data: { source: h5ad, path: data/sample.h5ad }

# d) a direct URL to a .h5ad or 10x .h5
data: { source: url,  url: https://example.org/dataset.h5ad }

# e) a dense CSV/TSV count matrix (optionally .gz), local or by URL
#    genes_are_rows: true (default) => rows are genes, columns are cells
data: { source: csv, url: https://.../counts.csv.gz, genes_are_rows: true }
```

GEO layouts vary; the loader tries `.h5ad` → 10x `.h5` → 10x mtx triplet (and
un-tars archives). If it can't auto-detect the matrix, download the supplementary
file manually and use `source: 10x_mtx | h5 | h5ad`, or `source: csv` if the
series ships a plain count table.

For **mouse** data, set `qc.mito_prefix: "mt-"`.

## Real public example: human endometrium (GSE111976)

A second worked example runs the whole pipeline on **real, public** single-cell
data — Wang et al. 2020, *Single-cell RNA-seq of human endometrium across the
natural menstrual cycle* ([GSE111976](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE111976)).
GEO distributes the counts as one gzipped CSV (~12 MB), loaded via `source: csv`:

```bash
python scrna_pipeline.py --config config/examples/endometrium_GSE111976.yaml
```

No access request is needed — it downloads straight from GEO. Outputs land in
`results/endometrium_GSE111976/` and `figures/endometrium_GSE111976/`. (If the
matrix loads transposed, set `genes_are_rows: false` in that config.)

## Starting from raw FASTQ (incl. controlled-access data)

This pipeline starts from a **count matrix**. Many datasets — including
controlled-access EGA studies (e.g. `EGAD*`) and SRA/ENA runs — are distributed
as raw **FASTQ**, which must be aligned and counted first:

1. **Obtain the data.** For controlled-access EGA datasets, submit a Data Access
   request to the study's DAC, sign the Data Access Agreement, then download with
   [`pyega3`](https://github.com/EGA-archive/pyega3):
   `pyega3 -cf credentials.json fetch EGAD50000001017`.
2. **Align + count** to produce a cell × gene matrix. For 10x data use
   **Cell Ranger** (`cellranger count`) or **STARsolo**; the output
   `filtered_feature_bc_matrix/` (or `.h5`) is what you feed in next.
3. **Run this pipeline** on the resulting matrix:
   `data: { source: 10x_mtx, path: <sample>/outs/filtered_feature_bc_matrix }`
   (or `source: h5`).

> **Keep controlled-access data private.** Data used under a DAA (e.g. EGA
> datasets) is personal data — never commit the raw data, matrices, or
> identifiable outputs to a public repo. Run these in a **private** repository
> and follow the agreement's terms.

## Parameters

All analysis choices live in `config/config.yaml` — QC thresholds (`n_mads`,
`min_genes`, `min_cells`, mito cap), `n_top_genes`, `n_pcs`, `leiden_resolution`,
marker `method`, and the `random_state`. Change them there, not in the code, so a
run is fully described by its config file.

## Reproducibility

- Pinned dependency versions (`environment.yml`, `requirements.txt`).
- Fixed random seeds for PCA, neighbors, UMAP, Leiden, and scrublet.
- One config file fully specifies a run; commit the config alongside results.

## Layout

```text
scrna_pipeline.py        the pipeline (load / QC / preprocess / cluster / markers)
config/config.yaml       dataset selection + all parameters
config/examples/         ready-to-run configs (e.g. endometrium_GSE111976.yaml)
notebooks/walkthrough.ipynb   interactive, cell-by-cell version
environment.yml          pinned conda environment
requirements.txt         pip alternative
Makefile                 make setup / make run / make clean
results/  figures/       outputs (git-ignored except .gitkeep)
```

## License

MIT — see `LICENSE`.
