# UK DRI scdownstream: Changelog

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/)
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## Unreleased — UK DRI fork

Forked from [nf-core/scdownstream](https://github.com/nf-core/scdownstream) at upstream `dev` commit
`fb5a421` (September 2025) and developed independently since. See
[`ACKNOWLEDGEMENTS.md`](ACKNOWLEDGEMENTS.md) for provenance. The full commit history is available via
`git log dev..dev_ukdri`.

### `Added`

- Three sequential entry points, replacing the single-pass workflow: `-entry qc_clustering`,
  `-entry downstream` and `-entry differential_genes`, chained by passing each stage's
  `<name>_qc_clustering.h5ad` / `<name>_downstream.h5ad` to the next as `--base_adata`.
  The single-pass workflow (`workflows/scdownstream.nf`) has been removed; running without
  `-entry` runs `qc_clustering` with a warning.
- Pseudobulk differential expression: decoupler pseudobulk aggregation, count/cell filtering, split
  per group label, and PyDESeq2 per group label × contrast. Contrasts are read from a TSV following
  the nf-core/differentialabundance definition.
- Quarto reporting, replacing the notebook-based reports: a QC/clustering report, a downstream
  analysis report and a differential expression report, with searchable tables capped by
  `--report_table_row_limit`.
- Automatic cell filtering from N-MAD outlier thresholds (`--automatic_cell_filtering`), plus
  explicit global QC filter parameters.
- Multi-resolution Leiden clustering in a single step (`--clustering_resolutions`), writing one
  `leiden_<res>` column per resolution.
- Gene set enrichment (`--enrich_*`) and marker gene export to JSON (`--markers_*`).
- Public container images for the report, enrichment and differential expression modules:
  `docker.io/nhecker/scanpy-report:1.11.4-coreinf0.4` and `docker.io/nhecker/pydeseq2:0.1`, built from Dockerfiles in
  the repository.
- `--memory_scale`, which scales every memory request in `conf/base.config`.
- `--qc_only`, to stop after per-sample QC and cell type annotation.
- Per-sample `n_hvgs` and `automatic_cell_filtering` columns in the samplesheet schema.
- `--metadata` / `--metadata_sample_col` for `-entry qc_clustering`: per-sample metadata from a TSV,
  added to `obs` before QC by the new `ADATA_ADDMETADATA` module. It stops on a missing sample,
  duplicate ids or a column that already exists in `obs`. The columns are kept through the merge,
  also with `--base_adata`.
- `--filtering_keep_outliers`: cells that fail the QC thresholds are kept and marked in a bool
  `obs["outlier"]` column instead of removed; the QC/clustering report shows them (UMAP,
  per-sample table, QC-metric violins).
- `--doublet_removal`: remove cells called doublets by at least `--doublet_detection_threshold`
  methods right after detection (default off: doublets are only annotated, as before).
  `--doublet_detection_threshold` is no longer inactive.
- `--pca_n_comps` (PCs computed, default 50) and `--neighbors_n_pcs` (PCs used for the PCA
  neighbour graph and the PCA UMAP, default all), plus an elbow plot in the QC/clustering report
  to choose them. With the defaults the results are unchanged.
- `--umap_color_by` / `--umap_color_by_embeddings`: both Quarto reports plot the named `obs`
  columns on each listed UMAP, by default the PCA and scVI UMAPs (before and after integration).
- `ADATA_ADDSAMPLE` sets the `sample` column in `obs` right after loading, instead of only in
  `ADATA_UNIFY`. A differing existing `sample` column is kept as `sample_original`.

### `Changed`

- Stage outputs are named after `--name` and the stage: `<name>_qc_clustering.h5ad`/`.rds`,
  `<name>_downstream.h5ad`/`.rds`, `<name>_downstream_markers.json.gz`, and the reports
  `<name>_qc_clustering_report.html` / `<name>_downstream_report.html` (previously
  `integrated_scvi_finalized.*` for stage 1 and `<name>_finalized.*` for stage 2).
- The curated tool set is now scrublet for doublet detection and scVI for integration; other tools
  remain in the codebase pending curation and validation.
- Doublet detection now runs before ambient RNA correction and filtering.
- LIANA+ uses local HCOP ortholog tables (`--ortholog_hcop_directory`, no default) for non-human data.
- `rank_genes_groups` uses the Wilcoxon test.
- `--integration_hvgs` default raised from 0 to 5000.
- `--unify_gene_symbols` is no longer supported: HUGO-based unification only applies to human data.
  Gene symbol resolution, duplicate handling and isoform aggregation are unaffected.
- Most intermediate outputs are now published only when `--save_intermediates` is set; the finalized
  objects are published at the top level of `--outdir`.
- Cell type predictions are merged into the per-sample objects as `obs` columns during finalisation.
- `ADATA_UNIFY` only copies an existing `sample` column to `sample_original` when its values differ
  from the sample id. Before, its check was always true, so any `sample` column was copied.
- `ADATA_MERGE` always keeps the `sample` column (and the `--metadata` columns). Before, with a
  `--base_adata` that had no `sample` column, it was dropped from the merged object.
- Documentation rewritten for the fork: README, `docs/usage.md`, `docs/output.md`, a new
  `ACKNOWLEDGEMENTS.md`, and an updated `CITATIONS.md`. Cluster-specific guidance lives on the
  UK DRI Informatics wiki rather than in the repository.

### `Known issues`

See [Status and known limitations](README.md#changes-and-known-limitations).

---

## v0.0.1dev - [2024-10-17]

Initial release of nf-core/scdownstream, created with the [nf-core](https://nf-co.re/) template.

### `Added`

- Added `singleR` module for automated cell type annotation.

### `Fixed`

### `Dependencies`

### `Deprecated`
