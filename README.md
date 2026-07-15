# scRNA-seq-reproducible

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
```

GEO layouts vary; the loader tries `.h5ad` → 10x `.h5` → 10x mtx triplet (and
un-tars archives). If it can't auto-detect the matrix, download the supplementary
file manually and use `source: 10x_mtx | h5 | h5ad`.

For **mouse** data, set `qc.mito_prefix: "mt-"`.

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
scrna_pipeline.py     the pipeline (load / QC / preprocess / cluster / markers)
config/config.yaml    dataset selection + all parameters
environment.yml       pinned conda environment
requirements.txt      pip alternative
Makefile              make setup / make run / make clean
results/  figures/    outputs (git-ignored except .gitkeep)
```

## License

MIT — see `LICENSE`.
