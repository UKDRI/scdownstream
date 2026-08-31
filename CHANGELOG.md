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
  `<name>_finalized.h5ad` to the next as `--base_adata`.
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
- `--singularity_cache_dir` for resolving locally built Apptainer/Singularity `.sif` images instead
  of pulling remote containers; defaults to `$NXF_SINGULARITY_CACHEDIR` / `$NXF_APPTAINER_CACHEDIR`.
- `--memory_scale`, which scales every memory request in `conf/base.config`.
- `--qc_only`, to stop after per-sample QC and cell type annotation.
- Per-sample `n_hvgs` and `automatic_cell_filtering` columns in the samplesheet schema.

### `Changed`

- The curated tool set is now scrublet for doublet detection and scVI for integration; other tools
  remain in the codebase pending curation and validation.
- Doublet detection now runs before ambient RNA correction and filtering.
- LIANA+ uses local HCOP ortholog tables (`--ortholog_hcop_directory`) for non-human data.
- `rank_genes_groups` uses the Wilcoxon test.
- `--integration_hvgs` default raised from 0 to 5000.
- `--unify_gene_symbols` is no longer supported: HUGO-based unification only applies to human data.
  Gene symbol resolution, duplicate handling and isoform aggregation are unaffected.
- Most intermediate outputs are now published only when `--save_intermediates` is set; the finalized
  objects are published at the top level of `--outdir`.
- Cell type predictions are merged into the per-sample objects as `obs` columns during finalisation.
- Documentation rewritten for the fork: README, `docs/usage.md`, `docs/output.md`, a new
  `ACKNOWLEDGEMENTS.md`, and an updated `CITATIONS.md`. Cluster-specific guidance lives on the
  UK DRI Informatics wiki rather than in the repository.

### `Known issues`

See [Status and known limitations](README.md#status-and-known-limitations).

---

## v0.0.1dev - [2024-10-17]

Initial release of nf-core/scdownstream, created with the [nf-core](https://nf-co.re/) template.

### `Added`

- Added `singleR` module for automated cell type annotation.

### `Fixed`

### `Dependencies`

### `Deprecated`
