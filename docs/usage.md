# UK DRI scdownstream: Usage

> This is a UK DRI fork of nf-core/scdownstream and is no longer part of the nf-core community.
> Parameter documentation is generated from [`nextflow_schema.json`](../nextflow_schema.json); run
> `nextflow run <pipeline> --help` to list every parameter and its default. There is no
> corresponding page on nf-co.re.

## Contents

- [Choosing a stage](#choosing-a-stage)
- [Samplesheet input](#samplesheet-input)
- [Filtered and unfiltered matrices](#filtered-and-unfiltered-matrices)
- [Running the pipeline](#running-the-pipeline)
- [`qc_clustering`](#qc_clustering)
- [`downstream`](#downstream)
- [`differential_genes`](#differential_genes)
- [Supported tool choices](#supported-tool-choices)
- [Container images](#container-images)
- [Reference data](#reference-data)
- [Cell type annotation](#cell-type-annotation)
- [Ambient RNA correction](#ambient-rna-correction)
- [Reference mapping](#reference-mapping)
- [GPU acceleration](#gpu-acceleration)
- [Core Nextflow arguments](#core-nextflow-arguments)
- [Custom configuration](#custom-configuration)

## Choosing a stage

The pipeline is split into three sequential stages, each selected with the `--stage` parameter.
**`--stage` is always required.** Each stage consumes the previous stage's finalized `.h5ad` via
`--base_adata`.

| Stage | `--stage`            | Required input                                          | Main output                                       |
| ----- | -------------------- | ------------------------------------------------------- | ------------------------------------------------- |
| 1     | `qc_clustering`      | `--input samplesheet.csv`                               | `<name>_qc_clustering.h5ad`                       |
| 2     | `downstream`         | `--base_adata` (stage 1 h5ad)                           | `<name>_downstream.h5ad`, `<name>_downstream_markers.json.gz` |
| 3     | `differential_genes` | `--base_adata` (stage 2 h5ad) + `--diffgenes_contrasts` | per-contrast DE tables                            |

`--name` sets the identifier used in output file names; it defaults to the stage name.

> [!IMPORTANT]
> Without a valid `--stage` the pipeline stops with an error before any task runs. Nextflow's
> `-entry` flag is no longer used: Nextflow 26.04's strict syntax parser does not support it. The
> upstream single-pass workflow has been removed.

## Samplesheet input

The samplesheet is required by `--stage qc_clustering`. It is a comma-separated file with a header row
and at least two columns: `sample`, and at least one of `filtered` / `unfiltered`.

```bash
--input '[path to samplesheet file]'
```

### Minimal samplesheet

```csv title="samplesheet.csv"
sample,unfiltered
sample1,/absolute/path/to/sample1.h5ad
sample2,relative/path/to/sample2.rds
sample3,/absolute/path/to/sample3.csv
```

### Full samplesheet

Optional columns enable more advanced features:

```csv title="samplesheet.csv"
sample,filtered,unfiltered,batch_col,label_col,unknown_label,expected_cells,ambient_correction,ambient_corrected_integration
sample1,/absolute/path/to/sample1_filtered.h5ad,/absolute/path/to/sample1.h5ad,batch,cell_type,unknown,5000,true,false
sample2,relative/path/to/sample2_filtered.rds,relative/path/to/sample2.rds,batch_id,annotation,unannotated,3000,false,
sample3,/absolute/path/to/sample3_filtered.csv,/absolute/path/to/sample3.csv,,,,,true,true
```

For CSV input files, `batch_col`, `label_col` and `unknown_label` have no effect, as no additional
metadata is available in a CSV.

| Column                          | Description                                                                                                                                                                                                                                                                                             |
| ------------------------------- | ------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- |
| `sample`                        | Unique sample identifier. Added to the output objects as a `sample` column.                                                                                                                                                                                                                             |
| `filtered`                      | Path to an `h5ad`, `h5`, `rds` or `csv` file. `rds` files may hold any object convertible to a `SingleCellExperiment` via [Seurat `as.SingleCellExperiment`](https://satijalab.org/seurat/reference/as.singlecellexperiment). `csv` files should hold a matrix with genes as columns and cells as rows. |
| `unfiltered`                    | As `filtered`, but for the unfiltered CellRanger or nf-core/scrnaseq output. If not provided, only `decontX` can be used for ambient RNA removal.                                                                                                                                                       |
| `batch_col`                     | Column in the input file holding batch information. If absent, the whole object is treated as one batch. Renamed to `batch` during execution.                                                                                                                                                           |
| `symbol_col`                    | Where this sample's gene symbols are. **Set per sample; there is no global `--symbol_col` parameter.** Empty or absent: `index` (the matrix row names / var index are the symbols). `none`: convert the gene IDs in `geneid_col` to symbols via MyGene.info. Any other value: the name of a var column holding the symbols. See [Gene symbols](#gene-symbols). |
| `geneid_col`                    | Var column holding gene IDs (e.g. Ensembl). Only used when `symbol_col` is `none`. Empty or absent: `index`.                                                                                                                                                                                         |
| `label_col`                     | Column holding cell type information. Defaults to `label`. If absent, the pipeline creates it and fills it with `unknown`. Renamed to `label` during execution.                                                                                                                                         |
| `unknown_label`                 | Value in `label_col` to treat as unknown. Defaults to `unknown`. Renamed to `unknown` during execution.                                                                                                                                                                                                 |
| `counts_layer`                  | Layer holding the raw counts matrix. Defaults to `X`.                                                                                                                                                                                                                                                   |
| `expected_cells`                | Expected number of cells, passed to CellBender for empty droplet detection.                                                                                                                                                                                                                             |
| `ambient_correction`            | Whether to run ambient RNA correction for this sample (`true` uses the globally configured method). Defaults to `true`.                                                                                                                                                                                 |
| `ambient_corrected_integration` | Whether to use ambient-corrected counts for integration for this sample, overriding the global `--ambient_corrected_integration`.                                                                                                                                                                       |
| `n_hvgs`                        | Number of highly variable genes for this sample. Defaults to `3000`.                                                                                                                                                                                                                                    |

### Gene symbols

The pipeline needs gene **symbols** (e.g. `CD3E`, `MT-CO1`), not only gene IDs: CellTypist models,
LIANA+ ligand–receptor resources and gene set enrichment all match on symbols. Where a sample's
symbols are is set **per sample** with the samplesheet column `symbol_col`. There is **no global
parameter** for it (`--symbol_col` was removed); a sample that leaves the column empty, or a
samplesheet without the column, uses `index`.

| Input                                                                  | `symbol_col`                | `geneid_col`  |
| ---------------------------------------------------------------------- | --------------------------- | ------------- |
| Symbols are the matrix row names / var index (e.g. CellRanger `.h5`)    | empty or `index`            | —             |
| Symbols are in a var column, e.g. `gene_symbols`                       | `gene_symbols`              | —             |
| Only gene IDs, e.g. Ensembl IDs in the index                           | `none`                      | empty or `index` |
| Only gene IDs, in a var column, e.g. `gene_ids`                        | `none`                      | `gene_ids`    |

```csv title="samplesheet.csv"
sample,filtered,symbol_col,geneid_col
sample1,/path/to/sample1_cellranger.h5,,
sample2,/path/to/sample2.h5ad,gene_symbols,
sample3,/path/to/sample3_ensembl_only.h5ad,none,gene_ids
```

> [!WARNING]
> `symbol_col = none` has three limits:
> - it looks the symbols up online at MyGene.info, so the compute nodes need internet access;
> - the lookup is for **human** genes only;
> - the conversion runs when the samples are unified, after per-sample QC and CellTypist. CellTypist
>   therefore cannot use it: with `--celltypist_model` set, a `none` sample stops with
>   "Symbol column none not found".
>
> For non-human data, or with CellTypist, provide input that already has a symbol column.

> [!IMPORTANT]
> Mitochondrial QC does **not** use `symbol_col`: it flags `MT-`/`mt-` genes in the var index, or
> in a var column named `gene_name`, `gene_names`, `gene_symbol` or `gene_symbols` if one exists.
> For input with Ensembl IDs in the index and none of these columns, the mitochondrial percentage is
> 0 for every cell and `max_mito_percentage` removes nothing.

### Per-sample QC thresholds

The cell and gene filtering thresholds can be set per sample. Each column overrides one global
parameter for that sample:

| Column                     | Overrides                    | Meaning                                                    |
| -------------------------- | ---------------------------- | ---------------------------------------------------------- |
| `min_genes`                | `--min_genes`                | Minimum genes per cell                                     |
| `max_genes`                | `--max_genes`                | Maximum genes per cell                                     |
| `min_counts_cell`          | `--min_counts`               | Minimum counts per cell                                    |
| `max_counts`               | `--max_counts`               | Maximum counts per cell                                    |
| `max_mito_percentage`      | `--max_mito_percentage`      | Maximum percentage of mitochondrial counts per cell        |
| `min_cells`                | `--min_cells`                | Minimum cells per gene (gene filtering)                    |
| `min_counts_gene`          | `--min_counts_gene`          | Minimum counts per gene (gene filtering)                   |
| `automatic_cell_filtering` | `--automatic_cell_filtering` | `true` / `false`: use N-MAD thresholds for this sample     |

How the threshold for a sample is chosen:

1. **A value in the samplesheet wins**, including `0` and `false`.
2. **An empty cell (or a missing column) falls back to the global parameter.** A samplesheet
   without these columns behaves exactly as before.
3. **With automatic filtering on** (per sample or globally), the N-MAD thresholds replace the fixed
   ones (`max_genes`, `max_counts`, `max_mito_percentage`), but `min_genes` and `min_counts_cell`
   (or `--min_genes` / `--min_counts`) still apply as lower floors. See
   [Automatic thresholds](#automatic-thresholds).
4. **Gene filtering** (`min_cells`, `min_counts_gene`) applies in both modes.

For example, to filter one low-quality sample more strictly and let another use automatic
thresholds while everything else uses the global settings:

```csv title="samplesheet.csv"
sample,filtered,min_genes,max_mito_percentage,automatic_cell_filtering
sample1,/path/to/sample1.h5ad,500,10,
sample2,/path/to/sample2.h5ad,,,true
sample3,/path/to/sample3.h5ad,,,
```

The thresholds applied to each sample are reported in its MultiQC section.

An [example samplesheet](../assets/samplesheet.csv) is provided with the pipeline.

## Filtered and unfiltered matrices

`unfiltered` matrices still contain empty droplets; `filtered` matrices have had empty droplets
removed. A more technical definition is
[here](https://support.10xgenomics.com/single-cell-gene-expression/software/pipelines/latest/output/matrices).
CellRanger provides both; most other quantification tools provide only the unfiltered matrix.

The pipeline handles three cases:

1. **Both matrices available** — provide both. The unfiltered matrix is used for ambient RNA removal
   and the filtered matrix for everything else.
2. **Only the filtered matrix** — provide it and it is used for all steps. Only `decontx` can be used
   for ambient RNA removal, as the other methods need the unfiltered matrix.
3. **Only the unfiltered matrix** — provide it and the pipeline creates a filtered matrix by
   identifying empty droplets with CellBender.

## Running the pipeline

```bash
nextflow run UKDRI/scdownstream -r dev_ukdri --stage qc_clustering \
   --input ./samplesheet.csv \
   --outdir ./results \
   -profile apptainer
```

Nextflow creates the following in your working directory:

```bash
work                # Nextflow working files
<OUTDIR>            # Results, at the location given by --outdir
.nextflow_log       # Nextflow log file
# plus other hidden Nextflow files, e.g. run history and old logs
```

For repeated runs, put the parameters in a YAML or JSON file and pass it with `-params-file`:

```bash
nextflow run UKDRI/scdownstream -r dev_ukdri --stage qc_clustering -profile apptainer -params-file params.yaml
```

```yaml title="params.yaml"
input: "./samplesheet.csv"
outdir: "./results/"
name: "my_study"
species: "human"
```

> [!WARNING]
> Do not use `-c <file>` to specify parameters — this will result in errors. Config files given with
> `-c` must only be used for process resource specifications, other infrastructural tweaks, or module
> arguments (`ext.args`).

One params file per stage — `params_qc_clustering.yml`, `params_downstream.yml`,
`params_differential_genes.yml` — keeps the three runs reproducible and easy to re-run
independently.

## `qc_clustering`

Stage 1: per-sample QC, cell type annotation, merging, scVI integration,
embeddings and Leiden clustering.

```bash
nextflow run UKDRI/scdownstream -r dev_ukdri --stage qc_clustering \
   -profile apptainer \
   --input samplesheet.csv \
   --name my_study \
   --species human \
   --outdir results/qc_clustering
```

Produces `results/qc_clustering/my_study_qc_clustering.h5ad` (plus `.rds`), a Quarto QC/clustering report
under `report/`, and a MultiQC report.

| Parameter                         | Default          | Description                                                                                                                                                            |
| --------------------------------- | ---------------- | ---------------------------------------------------------------------------------------------------------------------------------------------------------------------- |
| `--input`                         | —                | Samplesheet (required).                                                                                                                                                |
| `--metadata`                      | `null`           | Tab-separated per-sample metadata. Every column except `--metadata_sample_col` is added to the cells' `obs` — see [Per-sample metadata](#per-sample-metadata).         |
| `--metadata_sample_col`           | `sample`         | Column of the `--metadata` file holding the samplesheet `sample` ids.                                                                                                  |
| `--umap_color_by`                 | `''`             | Comma-separated `obs` columns to plot on the UMAPs in the report — see [UMAP colourings](#umap-colourings-in-the-reports). |
| `--umap_color_by_embeddings`      | `X_umap_pca,X_umap_scvi`| UMAPs those columns are plotted on: before (PCA) and after (scVI) integration. |
| `--name`                          | `qc_clustering`  | Identifier used in output file names.                                                                                                                                  |
| `--species`                       | `human`          | `human` or `mouse`. Set explicitly for mouse data.                                                                                                                     |
| `--ambient_correction`            | `none`           | `none`, `decontx`, `cellbender`, `soupx` or `scar`. Off by default, because the expected input is already CellBender-corrected. See [Ambient RNA correction](#ambient-rna-correction). |
| `--ambient_corrected_integration` | `false`          | Use ambient-corrected counts for integration rather than storing them as extra layers.                                                                                 |
| `--doublet_detection`             | `scrublet`       | Doublet detection method, or `none` to skip. See [Supported tool choices](#supported-tool-choices).                                                                    |
| `--doublet_removal`               | `false`          | Remove cells called doublets by at least `--doublet_detection_threshold` methods, right after doublet detection. Off: doublets are only annotated. See [Doublets and QC outliers](#doublets-and-qc-outliers). |
| `--doublet_detection_threshold`   | `1`              | With `--doublet_removal`: number of `--doublet_detection` methods that must call a cell a doublet. |
| `--automatic_cell_filtering`      | `false`          | Derive per-sample thresholds from N-MAD outlier detection instead of the fixed thresholds below. `--min_genes` / `--min_counts` still apply as lower floors; consider setting them very low in this mode. See [Automatic thresholds](#automatic-thresholds). Can be set per sample in the samplesheet. |
| `--filtering_keep_outliers`       | `false`          | Keep cells that fail the QC thresholds and mark them in `obs["outlier"]` instead of removing them. Cells left with no counts after gene filtering are removed regardless. See [Doublets and QC outliers](#doublets-and-qc-outliers). |
| `--min_genes`                     | `200`            | Minimum genes per cell. Also the floor for the automatic lower bound with `--automatic_cell_filtering`. Can be set per sample in the samplesheet. |
| `--max_genes`                     | `false`          | Maximum genes per cell (`false` disables). Can be set per sample in the samplesheet. |
| `--min_cells`                     | `5`              | Minimum cells per gene. Can be set per sample in the samplesheet. |
| `--min_counts`                    | `0`              | Minimum counts per cell. Also the floor for the automatic lower bound with `--automatic_cell_filtering`. Can be set per sample in the samplesheet. |
| `--max_counts`                    | `false`          | Maximum counts per cell (`false` disables). Can be set per sample in the samplesheet. |
| `--min_counts_gene`               | `0`              | Minimum counts per gene. Can be set per sample in the samplesheet. |
| `--max_mito_percentage`           | `25`             | Maximum percentage of mitochondrial counts per cell. Can be set per sample in the samplesheet. |
| `--duplicate_var_resolution`      | `sum`            | How to resolve duplicate gene symbols: `mean`, `sum`, `max` or `make_unique`.                                                                                          |
| `--aggregate_isoforms`            | `false`          | Aggregate isoform-level features.                                                                                                                                      |
| `--n_hvgs`                        | `3000`           | Highly variable genes used for the PCA/UMAP embedding.                                                                                                                 |
| `--pca_n_comps`                   | `50`             | Principal components computed and stored in `X_pca`. |
| `--neighbors_n_pcs`               | all              | Principal components used as input for the PCA neighbour graph (and the PCA UMAP), the first N of `X_pca`. See [Choosing the number of PCs](#choosing-the-number-of-pcs). |
| `--integration_methods`           | `scvi`           | Integration method. Keep `scvi` — see [Supported tool choices](#supported-tool-choices).                                                                               |
| `--integration_hvgs`              | `5000`           | Highly variable genes used for integration.                                                                                                                            |
| `--clustering_resolutions`        | `0.5,1.0`        | Comma-separated Leiden resolutions. One `leiden_<res>` column is written per resolution, and the **first** resolution is copied to `leiden` as the default clustering. |
| `--cluster_neighbors`             | `neighbors_scvi` | Neighbour graph used for clustering. The default assumes scVI integration.                                                                                             |
| `--celltypist_model`              | `null`           | Comma-separated [CellTypist](https://www.celltypist.org/models) model names, or a path to a custom `.pkl` model.                                                       |
| `--scvi_model` / `--scanvi_model` | `null`           | Pre-trained model for reference mapping — see [Reference mapping](#reference-mapping).                                                                                 |
| `--save_intermediates`            | `false`          | Publish intermediate objects as well as the final ones.                                                                                                                |

scVI itself is tuned with `--scvi_n_latent` (30), `--scvi_n_hidden` (128), `--scvi_n_layers` (2),
`--scvi_dispersion` (`gene`), `--scvi_gene_likelihood` (`zinb`), `--scvi_max_epochs`,
`--scvi_categorical_covariates` and `--scvi_continuous_covariates`.

> [!NOTE]
> `--species` defaults to `human`. Set it explicitly for mouse data: the human default is applied
> silently, and gene set enrichment and LIANA+ ortholog mapping would then use the wrong species.

### Doublets and QC outliers

Per-sample QC runs in this order: QC metrics, doublet detection, ambient RNA correction, then cell
filtering on the QC thresholds (fixed, or automatic N-MAD thresholds with
`--automatic_cell_filtering`). Two options decide what happens to the cells that are caught:

- **Doublets**: by default they are only annotated. scrublet writes `doublet_score` and
  `predicted_doublet` to `obs`, and no cells are removed. With `--doublet_removal`, cells called
  doublets by at least `--doublet_detection_threshold` of the `--doublet_detection` methods
  (default 1, i.e. any method) are removed right after detection, and the MultiQC cell-count table
  gets a "dedoubleted" row. Cells scrublet cannot score count as not doublets.
- **QC outliers**: by default, cells that fail the QC thresholds are removed. With
  `--filtering_keep_outliers` they are kept and marked in a bool `obs["outlier"]` column (`True` =
  fails at least one threshold). They then go through integration, clustering and the reports, so
  you can see where they fall on the UMAP and in which clusters before deciding to remove them. The
  QC/clustering report then shows a "QC outliers" section: a UMAP coloured by `outlier`, a table of
  outliers per sample and violin plots of the QC metrics. Gene filtering (`--min_cells`,
  `--min_counts_gene`) still applies, counted over all kept cells, outliers included.
- **Empty cells**: gene filtering (`--min_cells`, `--min_counts_gene`) runs after cell filtering, so a
  cell whose only genes are removed is left with no counts. Such cells are **always removed**, also
  with `--filtering_keep_outliers`: an empty barcode cannot be normalised or integrated and makes
  later steps fail. The number removed is in the MultiQC caption of each sample and in
  `uns["cell_filtering"]["n_empty_after_gene_filtering"]`.

The two options are independent: you can, for example, remove doublets and keep QC outliers.

#### Automatic thresholds

With `--automatic_cell_filtering`, each sample gets its own thresholds from the median and the median
absolute deviation (MAD) of its QC metrics:

| Metric                               | Lower bound                                     | Upper bound      |
| ------------------------------------ | ----------------------------------------------- | ---------------- |
| genes per cell (`n_genes_by_counts`) | median − 2.5 × MAD, but at least `--min_genes`  | median + 5 × MAD |
| counts per cell (`total_counts`)     | median − 2.5 × MAD, but at least `--min_counts` | median + 5 × MAD |
| mitochondrial % (`pct_counts_mt`)    | none                                            | median + 5 × MAD |

`--min_genes` and `--min_counts` (or a sample's `min_genes` / `min_counts_cell` column) are
**floors**: the automatic lower bound never goes below them. Automatic filtering can also be
switched on or off per sample with the `automatic_cell_filtering` column; see
[Per-sample QC thresholds](#per-sample-qc-thresholds). On
samples with a wide or skewed distribution, median − 2.5 × MAD is often below zero, and without a
floor no cell would be removed at the low end, near-empty barcodes included. With the defaults, cells
need at least 200 genes in automatic mode too. The fixed `--max_genes`, `--max_counts` and
`--max_mito_percentage` are not used in automatic mode.

> [!TIP]
> When you use `--automatic_cell_filtering`, consider setting `--min_genes` and `--min_counts` very
> low (e.g. `--min_genes 10 --min_counts 20`). The N-MAD thresholds then decide where the lower
> cut-off lies for each sample, and the floors only stop empty and near-empty barcodes from passing
> when the N-MAD bound drops to zero. With the default `--min_genes 200`, the floor often overrides
> the automatic bound, and samples with many small cells (e.g. nuclei with low RNA content) lose
> more cells than the N-MAD rule alone would remove.

The thresholds that were applied are in the per-sample MultiQC section (which also says when a floor
replaced the N-MAD value) and in `uns["n_mads_obs_thresholds"]`: `thresholds` holds the applied
values, `n_mads_thresholds` the raw N-MAD values and `floors` the floors.

> [!NOTE]
> Cells with no counts in the highly variable genes used for integration are also dropped before
> scVI, also with `--filtering_keep_outliers`, so a few cells with very few genes can be missing
> from the merged object.

To drop the marked outliers later, subset the object, e.g. `adata = adata[~adata.obs["outlier"]].copy()`.

### Per-sample metadata

Right after loading, before QC, every sample's cells get a `sample` column in `obs` holding the
samplesheet `sample` id. If the input object already has a `sample` column with other values, they
are kept as `sample_original` and a warning is printed; if `sample_original` also exists, the
pipeline stops rather than overwrite it.

To add further per-sample annotations (donor, diagnosis, sex, age, ...), pass `--metadata` a
tab-separated file with a header row and one row per sample:

```tsv
sample	donor	diagnosis	sex	age
sample1	D01	AD	F	71
sample2	D02	control	M	68
```

- The column named by `--metadata_sample_col` (default `sample`) holds the samplesheet `sample`
  ids. It is matched against the `sample` column in `obs`, and is not added itself.
- Every other column is added to `obs` for all cells of that sample. Numeric columns stay numeric;
  everything else, including `true`/`false`, is stored as a categorical of strings. Empty cells
  become missing values.
- Rows for samples that are not in the run are ignored.

The pipeline stops with an error when:

- a sample in the run has no row in the file;
- a sample id is empty or appears twice;
- the header has duplicate or empty column names (often a trailing tab);
- a metadata column already exists in the input object's `obs`;
- a metadata column is named `batch`, `label`, `sample_original` or `outlier`. The pipeline
  creates these itself, so rename them, e.g. `batch` to `seq_batch`.

The columns are added before QC, so they are in every object from then on: the per-sample QC
objects, the merged and integrated object, and the stage outputs used by
`downstream` and `differential_genes`. That means you can use them as scVI covariates
(`--scvi_categorical_covariates`, `--scvi_continuous_covariates`), and as `variable` or `blocking`
in the [contrasts file](#the-contrasts-file). With `--base_adata`, the metadata columns are kept
even if the base object lacks them; its cells then get `unknown` (or `NaN` for numeric columns).

Also avoid names that later steps write to `obs`: the QC metrics (`n_genes_by_counts`,
`total_counts`, `pct_counts_mt`, ...), `doublet_score` / `predicted_doublet`,
`celltypist:<model>` and `leiden_<resolution>`. Those steps would replace your column without a
warning. Only `--stage qc_clustering` reads `--metadata`.

### Choosing the number of PCs

Two parameters control the PCA-based embedding, and they are different things:

- `--pca_n_comps` (default `50`) is how many principal components are **computed** and stored
  in `X_pca`. On a small object the pipeline computes fewer, at most one less than the number of
  cells or of genes used for PCA (the highly variable genes), and prints a warning.
- `--neighbors_n_pcs` (default: all computed PCs) is how many of those PCs are **used** to build
  the PCA neighbour graph, i.e. the first N columns of `X_pca`. The PCA UMAP (`X_umap_pca`) is
  computed from this graph. It cannot be larger than `--pca_n_comps`.

To pick a value, look at the elbow plot in the PCA section of the QC/clustering report. It shows
the variance explained per PC, with a line at the number of PCs the graph used. Then rerun stage 1
with the new value and `-resume`:

```bash
nextflow run UKDRI/scdownstream -r dev_ukdri --stage qc_clustering -resume \
   <same options as before> \
   --neighbors_n_pcs 20
```

With `-resume`, QC and scVI integration come from the cache; only the steps from the PCA graph
onwards are rerun. Add `--pca_n_comps` only if you want more than 50 PCs. That also reruns the PCA
itself.

These parameters change the PCA graph and the PCA UMAP only. The scVI graph and the scVI UMAP,
and the default clustering built on the scVI graph, stay the same unless you cluster on the PCA
graph (`--cluster_neighbors neighbors_pca`).

### UMAP colourings in the reports

`--umap_color_by` names `obs` columns to plot on the UMAPs, for example columns added with
`--metadata`:

```bash
--metadata sample_metadata.tsv --umap_color_by diagnosis,sex,age
```

Both reports (the QC/clustering report of `qc_clustering` and the analysis report of `downstream`)
then get a "Requested UMAP plots" section with one block per UMAP in
`--umap_color_by_embeddings`, each showing every requested column. The default,
`X_umap_pca,X_umap_scvi`, shows the columns **before integration** (PCA UMAP) and **after
integration** (scVI UMAP). This helps you judge what integration did: a technical variable such as
a sequencing batch should mix after integration, while a biological one such as a cell type
should still separate.

- Numeric columns get a colour scale and categorical ones a legend. Columns with more than 50
  categories are plotted without a legend.
- Columns or UMAPs that are not in the object are skipped, with a note in the report listing what
  is available. The report never fails because of them.
- `--umap_for_plots` is separate: it still picks the UMAP for all other plots in the reports.

## `downstream`

Stage 2: marker genes per cluster, gene set enrichment, and LIANA+ cell–cell communication.

```bash
nextflow run UKDRI/scdownstream -r dev_ukdri --stage downstream \
   -profile apptainer \
   --base_adata results/qc_clustering/my_study_qc_clustering.h5ad \
   --name my_study \
   --species human \
   --outdir results/downstream
```

> [!NOTE]
> `--input` is not needed here: this stage works entirely from `--base_adata`. Only `--outdir` is
> strictly required by the parameter schema, so check your other parameters carefully — a typo in a
> parameter name is caught, but an omitted one falls back to its default.

| Parameter                         | Default                 | Description                                                                                                                                             |
| --------------------------------- | ----------------------- | ------------------------------------------------------------------------------------------------------------------------------------------------------- |
| `--base_adata`                    | —                       | The finalized h5ad from `qc_clustering` (required).                                                                                                     |
| `--name`                          | `downstream`            | Identifier used in output file names.                                                                                                                   |
| `--selected_clustering`           | `leiden`                | The `obs` column to group cells by for marker genes, enrichment and LIANA+. Use e.g. `leiden_1.0` to pick a specific resolution, or a cell type column. |
| `--species`                       | `human`                 | `human` or `mouse`. Drives enrichment and LIANA+ ortholog mapping.                                                                                      |
| `--enrich_min_in_group_fraction`  | `0.25`                  | Minimum fraction of cells in the group expressing a gene for it to enter enrichment.                                                                    |
| `--enrich_min_fold_change`        | `1.0`                   | Minimum fold change for a gene to enter enrichment.                                                                                                     |
| `--enrich_max_out_group_fraction` | `0.5`                   | Maximum fraction of cells outside the group expressing a gene.                                                                                          |
| `--ortholog_hcop_directory`       | `null`                  | Directory of HCOP ortholog tables; required for LIANA+ on non-human data — see [Reference data](#reference-data).                                       |
| `--markers_uns_key`               | `rank_genes_groups`     | `uns` key holding the marker results to export.                                                                                                         |
| `--markers_thr_adj_pvalue`        | `0.05`                  | Adjusted p-value threshold for the exported markers.                                                                                                    |
| `--markers_n_top`                 | `100`                   | Number of top markers per group to export.                                                                                                              |
| `--markers_pct_nz`                | `0.1`                   | Minimum fraction of non-zero expressing cells for an exported marker.                                                                                   |
| `--markers_min_logfc`             | `0`                     | Minimum log fold change for an exported marker.                                                                                                         |
| `--markers_sep`                   | `--`                    | Separator between the parts of the set ids in both marker-gene JSONs (`genes--cluster--<cluster>`, `genes--<contrast>--<group>--<up\|down>`). `--` rather than `__`, because contrast ids often contain `__`. Applies to `downstream` and `differential_genes`; to change it, set it in a `-params-file`. |
| `--report_table_row_limit`        | `250`                   | Maximum rows shown per table in the HTML report.                                                                                                        |
| `--report_max_celltype_labels`    | `50`                    | Maximum number of CellTypist labels shown per model in the "Predicted labels" sections of both reports. Labels are ranked by the number of predicted cells; the UMAP colours the top N and draws the rest in grey as `other`, and the count table lists the top N plus an `other` row. Raise it to see more labels, at the cost of a longer legend. |
| `--umap_for_plots`                | `X_umap_scvi`           | UMAP used for the report plots: `X_umap_scvi` (scVI graph) or `X_umap_pca` (PCA graph). Also applies to the stage 1 report.                            |
| `--umap_color_by`                 | `''`                    | Comma-separated `obs` columns to plot on the UMAPs in both reports, e.g. `sex,diagnosis`. See [UMAP colourings](#umap-colourings-in-the-reports). |
| `--umap_color_by_embeddings`      | `X_umap_pca,X_umap_scvi`| UMAPs the `--umap_color_by` columns are plotted on, one block each: before (PCA) and after (scVI) integration. |

Marker genes are computed with `scanpy.tl.rank_genes_groups` using the **Wilcoxon** test.

## `differential_genes`

Stage 3: pseudobulk aggregation with decoupler, then PyDESeq2 differential expression for every
combination of group label and contrast.

```bash
nextflow run UKDRI/scdownstream -r dev_ukdri --stage differential_genes \
   -profile apptainer \
   --base_adata results/downstream/my_study_downstream.h5ad \
   --diffgenes_contrasts contrasts.tsv \
   --diffgenes_group_col cell_type \
   --diffgenes_sample_col sample \
   --name my_study \
   --outdir results/differential_genes
```

> [!IMPORTANT]
> The input h5ad must contain a **`counts` layer** holding raw counts. Pseudobulk aggregation sums
> that layer; normalised data cannot be used.

| Parameter                  | Default | Description                                                                                                                         |
| -------------------------- | ------- | ----------------------------------------------------------------------------------------------------------------------------------- |
| `--base_adata`             | —       | Input h5ad with a `counts` layer (required).                                                                                        |
| `--diffgenes_contrasts`    | —       | Contrasts TSV (required) — see below.                                                                                               |
| `--diffgenes_group_col`    | `''`    | `obs` column holding the group label to test within, e.g. a cluster or cell type column. One DE analysis is run per distinct value. |
| `--diffgenes_sample_col`   | `''`    | `obs` column identifying the biological sample. Pseudobulk profiles are formed per sample.                                          |
| `--diffgenes_min_counts`   | `1000`  | Drop pseudobulk samples with fewer total counts.                                                                                    |
| `--diffgenes_min_cells`    | `10`    | Drop pseudobulk samples aggregated from fewer cells.                                                                                |
| `--diffgenes_min_samples`  | `2`     | Minimum pseudobulk samples required **on each side** of a contrast.                                                                 |
| `--diffgenes_markers_n_top`| `100`   | Maximum genes per group label, contrast and direction in `<name>_differential_genes.json.gz`. |
| `--diffgenes_markers_thr_adj_pvalue`| `0.05`  | Adjusted p-value threshold for the marker-gene JSON. |
| `--diffgenes_markers_min_logfc`| `0`     | Minimum absolute log2 fold change for the marker-gene JSON. |
| `--markers_sep`            | `--`    | Separator in the set ids of the marker-gene JSON, e.g. `genes--AD_vs_ctrl--Microglia--up`. Shared with `downstream`. |
| `--report_table_row_limit` | `250`   | Maximum rows shown per table in the HTML report.                                                                                    |

### The contrasts file

A **tab-separated** file following the
[nf-core/differentialabundance](https://nf-co.re/differentialabundance) contrasts definition,
validated against [`assets/schema_contrasts.json`](../assets/schema_contrasts.json). An example is
provided at [`assets/contrasts.tsv`](../assets/contrasts.tsv).

| Column                   | Required | Description                                                                                                                     |
| ------------------------ | -------- | ------------------------------------------------------------------------------------------------------------------------------- |
| `id`                     | yes      | Contrast name, used in output file names.                                                                                       |
| `variable`               | yes      | `obs` column holding the condition to compare.                                                                                  |
| `reference`              | yes      | Value of `variable` treated as the reference (baseline) group.                                                                  |
| `target`                 | yes      | Value of `variable` treated as the target group. Log fold changes are target vs reference.                                      |
| `blocking`               | no       | Additional `obs` columns to include in the design as covariates. **Comma-separated**, because the file itself is tab-separated. |
| `exclude_samples_col`    | no       | Sample-level `obs` column used to drop samples from **this contrast only**, e.g. `donor_id`. See [Excluding samples](#excluding-samples-from-a-contrast). |
| `exclude_samples_values` | no       | Values of `exclude_samples_col` whose samples are dropped. **Comma-separated.**                                                 |

The design formula is `~ variable [+ blocking…]`, fitted with `refit_cooks=True`.

#### Excluding samples from a contrast

Set `exclude_samples_col` and `exclude_samples_values` to leave samples out of one contrast, for
example outlier donors, while the other contrasts still use them. Every pseudobulk sample whose
value in `exclude_samples_col` is one of the `exclude_samples_values` is dropped before the model
is fitted and before the `--diffgenes_min_samples` check:

```tsv title="contrasts.tsv"
id	variable	reference	target	blocking	exclude_samples_col	exclude_samples_values
AD_vs_ctrl	condition	ctrl	AD		donor_id	D07,D12
```

- The column must be **sample-level** metadata, i.e. constant within each sample (donor, batch,
  sex, …). Pseudobulk aggregation keeps only columns that are constant within each profile, so a
  column that varies between cells of a sample is not available and the run stops with an error
  naming the column. Columns added with `--metadata` are sample-level and work.
- The task log reports how many samples were excluded and kept. Values that match no sample are
  reported as a warning (check for typos), as is a column without values or values without a
  column; in those cases nothing is excluded.
- If the exclusion leaves fewer than `--diffgenes_min_samples` samples on either side, the
  contrast is skipped for that group label, and the skip message says that samples were excluded.

> [!WARNING]
> A contrast is **silently skipped** for a given group label if either side has fewer than
> `--diffgenes_min_samples` surviving pseudobulk samples. No error is raised and no output file is
> written for that combination. Check the DE report's coverage against the contrasts you expected,
> and consult `de_manifest.tsv` in the report inputs to see which combinations actually produced
> results.

Every column named in the contrasts file — `variable` and each `blocking` column — is automatically
carried through pseudobulk aggregation alongside `--diffgenes_group_col`, so those columns must exist
in the input object's `obs`.

## Supported tool choices

This fork focuses on a curated set of tools — the approaches we have validated on UK DRI data. Please
use the values below. Further tools will be enabled as they are curated and validated; the parameter
schema is broader than this list, and values outside it are not yet supported.

| Step                   | Parameter                                   | Supported values                                           |
| ---------------------- | ------------------------------------------- | ---------------------------------------------------------- |
| Doublet detection      | `--doublet_detection`                       | `scrublet` (default), `none`                               |
| Ambient RNA correction | `--ambient_correction`                      | `none` (default), `decontx`, `soupx`, `cellbender`, `scar` |
| Integration            | `--integration_methods`                     | `scvi` (default)                                           |
| Clustering             | `--clustering_resolutions`                  | Leiden, one or more resolutions                            |
| Cell type annotation   | `--celltypist_model`                        | CellTypist                                                 |

Two consequences of the current tool set are worth knowing before you interpret results:

- **Doublets are annotated, not removed, by default.** scrublet writes its scores and its
  `predicted_doublet` call into the object. Set `--doublet_removal` to remove them; see
  [Doublets and QC outliers](#doublets-and-qc-outliers).
- **scVI drives everything after the merge.** The embeddings, clustering and reports are built on the
  scVI latent space, and `--cluster_neighbors` defaults to `neighbors_scvi` accordingly. Keep `scvi`
  in `--integration_methods`; if it is omitted, the post-merge steps produce no output and no error.

**Inactive parameters.** The following are accepted but currently have no effect:
`--skip_liana`, `--skip_rankgenesgroups`,
`--pseudobulk`, `--pseudobulk_groupby_labels`, `--pseudobulk_min_num_cells`, `--cluster_per_label`,
`--cluster_global`, `--base_embeddings`, `--base_label_col`. `--unify_gene_symbols` is no longer
supported and should be left at its default —
HUGO-based gene symbol unification only applies to human data and is not reliable enough to
recommend.

Gene symbols are still harmonised across samples without it: `symbol_col` / `geneid_col` resolution
(including MyGene.info conversion when `symbol_col` is `none`), duplicate resolution via
`--duplicate_var_resolution`, and optional isoform aggregation via `--aggregate_isoforms` all run
regardless.

## Container images

Every process runs from a public container image, so nothing has to be built by hand. Nextflow
pulls each image on first use. Under the `singularity` or `apptainer` profile it converts the image
to a `.sif` in `$NXF_SINGULARITY_CACHEDIR`. Two images are specific to this fork. They are hosted on
Docker Hub and pulled the same way; their Dockerfiles are kept in the repository only as the recipe
they were built from:

| Image (Docker Hub)                        | Used by                                                                                            | Recipe                                                 |
| ----------------------------------------- | -------------------------------------------------------------------------------------------------- | ------------------------------------------------------ |
| `docker.io/nhecker/scanpy-report:1.11.4-coreinf0.4` | `SCANPY_GENERATE_REPORT`, `SCANPY_GENERATE_REPORT_QC`, `PYDESEQ2_GENERATE_REPORT`, `SCANPY_ENRICH` | `modules/local/scanpy/report/Dockerfile`               |
| `docker.io/nhecker/pydeseq2:0.1`                    | `DIFFERENTIAL_GENES_PER_CONTRAST`                                                                  | `modules/local/pydeseq2/differential_genes/Dockerfile` |

Both extend the public `gcfntnu/scanpy:1.11.4` image. `DECOUPLER_PSEUDOBULK`, `FILTER_PSEUDOBULK` and
`SCANPY_EXPORT_MARKERS` use that public image directly. The images are published for `linux/amd64`
only.

Prefer `-profile apptainer` over `-profile docker`. Docker Hub limits anonymous pulls, and Docker
pulls on every host, while Apptainer pulls each image once into `$NXF_SINGULARITY_CACHEDIR`. See the
note in the [README](../README.md#container-images).

## Reference data

**HCOP orthologs (LIANA+).** LIANA+ maps its human-derived ligand–receptor resource onto
non-human data using HCOP ortholog tables. For any `--species` other than human, download the
fifteen-column table for your species from the
[HGNC HCOP downloads](https://www.genenames.org/download/hcop/tsv/), e.g.

```bash
mkdir -p hcop
curl -o hcop/human_mouse_hcop_fifteen_column.txt.gz \
    https://storage.googleapis.com/public-download-files/hcop/human_mouse_hcop_fifteen_column.txt.gz
```

and pass the directory with `--ortholog_hcop_directory hcop`. LIANA+ reads
`<directory>/human_<species>_hcop_fifteen_column.txt.gz`, so keep the original file name.
`--ortholog_hcop_directory` has no default. `--stage downstream` stops with an error if it is
missing for a non-human species, and human data does not need it.

**CellTypist models** are downloaded at runtime unless `--celltypist_model` is given a local `.pkl`
path.

## Cell type annotation

### CellTypist

Specify one or more models with `--celltypist_model` (comma-separated). Available models are listed
at [celltypist.org](https://www.celltypist.org/models). A path containing `/` is treated as a custom
model file.

> [!NOTE]
> Cell type predictions are merged into the per-sample objects as `obs` columns during the
> finalisation step of `qc_clustering`.

## Ambient RNA correction

Ambient RNA correction removes contaminating RNA from cell-free droplets. It is **off by default**
(`--ambient_correction none`), because the expected input, nf-core/scrnaseq output, has usually
already been corrected with CellBender, and correcting it again adds nothing. For input that was not
corrected upstream (e.g. plain CellRanger output), select a method with `--ambient_correction`:
`decontx`, `cellbender`, `soupx` or `scar`.

```bash
nextflow run UKDRI/scdownstream -r dev_ukdri --stage qc_clustering \
   --ambient_correction cellbender --input samplesheet.csv --outdir results
```

With a method set, correction can be switched off for individual samples from the samplesheet
(with the default `none`, the column has no effect):

```csv title="samplesheet.csv"
sample,filtered,unfiltered,ambient_correction
sample1,/path/to/sample1_filtered.h5ad,/path/to/sample1.h5ad,true
sample2,/path/to/sample2_filtered.h5ad,/path/to/sample2.h5ad,false
```

By default the corrected counts are stored as an additional layer (e.g.
`ambient_corrected_decontx`) while the raw counts stay in `X`, so integration uses the raw counts and
the corrected counts are available for inspection. To use the corrected counts for integration
instead, set `--ambient_corrected_integration true` globally, or the
`ambient_corrected_integration` column per sample.

> [!WARNING]
> When `ambient_corrected_integration` is enabled, the corrected counts **replace** the raw counts in
> `X` and the original raw counts are no longer available.

## Reference mapping

`qc_clustering` supports mapping new samples into the latent space of an existing scVI or scANVI
model, via `--scvi_model` and `--scanvi_model`. With an scANVI model this also transfers cell type
annotations to the new samples.

- **Pre-trained scVI model** — pass its path to `--scvi_model` and keep `scvi` in
  `--integration_methods`. Only the genes present in the reference model are used, rather than
  intersecting highly variable genes.
- **Adding samples to an existing object** — additionally pass the previous run's h5ad to
  `--base_adata`; the new samples are aggregated onto it after integration.

> [!NOTE]
> `--scanvi_model` exists for label transfer from a pre-trained scANVI model, but scANVI is not yet
> part of the curated tool set (see [Supported tool choices](#supported-tool-choices)): it produces
> embeddings without driving the rest of the workflow.

## GPU acceleration

:::warning{title="Experimental feature"}
This is an experimental feature and may produce errors. Please report issues on the
[UKDRI/scdownstream issue tracker](https://github.com/UKDRI/scdownstream/issues).
:::

:::info{title="Prerequisites"}

- Tested with Docker, Singularity and Apptainer. Other container technologies may work but are
  untested. Conda is not supported.
- CUDA 12.0 or later is required.
- GPUs must have a
  [Compute Capability](https://docs.nvidia.com/cuda/cuda-c-programming-guide/index.html#compute-capabilities)
  of 7.0 or higher.

:::

Tools with GPU support:

- CellBender
- scvi-tools — scVI/scANVI, scAR, solo

Add the `gpu` profile to use CUDA-enabled environments. All GPU-capable processes carry the
`process_gpu` label. You must also ensure those tasks land on GPU nodes. On a Slurm cluster with a
dedicated GPU queue:

```bash
process {
  withLabel:process_gpu {
    queue = '<gpu-queue>'
    clusterOptions = '--gpus 1'
  }
}
```

:::tip
See the [Nextflow Slurm documentation](https://www.nextflow.io/docs/latest/executor.html#slurm).
Depending on your cluster you may need `--gpus-per-node=1` or `--gres=gpu:1` instead of `--gpus 1`.
:::

:::tip
If jobs land on the right nodes but the GPU is not used, you may need:
`singularity.runOptions = '--no-mount tmp --writable-tmpfs --nv --env CUDA_VISIBLE_DEVICES=$CUDA_VISIBLE_DEVICES --env ROCR_VISIBLE_DEVICES=$ROCR_VISIBLE_DEVICES --env ZE_AFFINITY_MASK=$ZE_AFFINITY_MASK --env NVIDIA_VISIBLE_DEVICES=$CUDA_VISIBLE_DEVICES'`

The first part (`--no-mount tmp --writable-tmpfs --nv`) is already set by the `gpu` profile.
:::

## Core Nextflow arguments

> [!NOTE]
> These are Nextflow options and use a _single_ hyphen. Pipeline parameters use a double hyphen.

### `-entry`

No longer used. The stage is selected with the `--stage` pipeline parameter; see
[Choosing a stage](#choosing-a-stage).

### `-profile`

Chooses a configuration profile. Profiles bundled with the pipeline select how software is provided:

- `apptainer` — [Apptainer](https://apptainer.org/); **recommended**, and the UK DRI default
- `singularity` — [Singularity](https://sylabs.io/docs/)
- `docker` — [Docker](https://docker.com/)
- `podman`, `shifter`, `charliecloud` — other container engines
- `conda` — [Conda](https://conda.io/docs/); a last resort, and not supported for GPU runs
- `gpu` — enables CUDA-capable environments, combine with a container profile
- `test` — a complete configuration for automated testing, including links to test data

Multiple profiles can be combined and are applied in order, so later profiles override earlier ones:
`-profile apptainer,gpu`.

> [!IMPORTANT]
> Use a container profile for reproducibility. If `-profile` is omitted the pipeline runs locally and
> expects every tool to be on `PATH`, which is not recommended.

> [!NOTE]
> `-profile test_offline` is no longer supported. Use `-profile test`.

The pipeline still dynamically loads institutional configurations from
[nf-core/configs](https://github.com/nf-core/configs). Since this fork is not an nf-core pipeline,
there is no pipeline-specific config there and no plan to add one — use `-c` with your own config
instead.

### `-resume`

Restarts a pipeline, reusing cached results for steps whose inputs are unchanged. Both file names and
contents must match. See
[this blog post](https://www.nextflow.io/blog/2019/demystifying-nextflow-resume.html). You can resume
a specific run with `-resume [run-name]`; `nextflow log` lists previous run names.

### `-c`

Specifies an additional config file. Use it for resources and infrastructure, never for parameters.

## Custom configuration

### Resource requests

Each process has default CPU, memory and time requests, defined in [`conf/base.config`](../conf/base.config)
and selected by process labels (`process_low`, `process_medium`, `process_high`, `process_gpu`, …).
Failed jobs are automatically retried with increased resources for a limited number of attempts.

To raise the memory for the heaviest processes, put this in a config file and pass it with `-c`:

```groovy title="custom.config"
process {
  withLabel:process_high {
    memory = 225.GB
  }
}
```

This fork also provides `--memory_scale`, which multiplies every memory request in
`conf/base.config` by the given factor — a quick way to scale the whole pipeline up for large
datasets.

### Custom containers

To override the container used by a single process, use a `withName` selector in a `-c` config:

```groovy title="custom.config"
process {
  withName:SCANPY_ENRICH {
    container = '/path/to/my_image.sif'
  }
}
```

The images this fork uses are listed under [Container images](#container-images).

### Custom tool arguments

Extra arguments can be passed to a process via `ext.args`:

```groovy title="custom.config"
process {
  withName:SCANPY_LEIDEN {
    ext.args = '--flavor igraph'
  }
}
```

### Updating the pipeline

Nextflow caches the pipeline code it pulls from GitHub. To make sure you are running the current
version:

```bash
nextflow pull UKDRI/scdownstream
```

### Reproducibility

Pin the revision you run with `-r`, e.g. `-r dev_ukdri` for the development branch or a specific
commit SHA or tag for a fixed version. The revision is recorded in the run reports. Sharing the
`-params-file` alongside it makes a run reproducible without a long command line.

> [!TIP]
> If you share a params file (for example as supplementary material), remove cluster-specific paths
> such as `--ortholog_hcop_directory`.

## Running in the background

Nextflow must keep running until the pipeline finishes. Use `-bg` to detach it from your terminal
(logs go to a file), or run it inside `screen`/`tmux`. On an HPC, the usual pattern is to run
Nextflow itself as a cluster job, from which it submits the individual tasks; note that its job
allocation sizes only the driver, not the work.

## Nextflow memory requirements

The Nextflow JVM can request a large amount of memory. Limiting it in your shell profile (typically
`~/.bashrc`) is recommended:

```bash
NXF_OPTS='-Xms1g -Xmx4g'
```
