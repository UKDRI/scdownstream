# UK DRI scdownstream: Usage

> This is a UK DRI fork of nf-core/scdownstream and is no longer part of the nf-core community.
> Parameter documentation is generated from [`nextflow_schema.json`](../nextflow_schema.json); run
> `nextflow run <pipeline> --help` to list every parameter and its default. There is no
> corresponding page on nf-co.re.

## Contents

- [Choosing an entry point](#choosing-an-entry-point)
- [Samplesheet input](#samplesheet-input)
- [Filtered and unfiltered matrices](#filtered-and-unfiltered-matrices)
- [Running the pipeline](#running-the-pipeline)
- [`qc_clustering`](#qc_clustering)
- [`downstream`](#downstream)
- [`differential_genes`](#differential_genes)
- [Supported tool choices](#supported-tool-choices)
- [Local container images](#local-container-images)
- [Reference data](#reference-data)
- [Cell type annotation](#cell-type-annotation)
- [Ambient RNA correction](#ambient-rna-correction)
- [Reference mapping](#reference-mapping)
- [GPU acceleration](#gpu-acceleration)
- [Core Nextflow arguments](#core-nextflow-arguments)
- [Custom configuration](#custom-configuration)

## Choosing an entry point

The pipeline is split into three sequential stages, each selected with Nextflow's `-entry` flag.
**`-entry` is always required.** Each stage consumes the previous stage's finalized `.h5ad` via
`--base_adata`.

| Stage | Entry point          | Required input                                          | Main output                                       |
| ----- | -------------------- | ------------------------------------------------------- | ------------------------------------------------- |
| 1     | `qc_clustering`      | `--input samplesheet.csv`                               | `<name>_finalized.h5ad`                           |
| 2     | `downstream`         | `--base_adata` (stage 1 h5ad)                           | `<name>_finalized.h5ad`, `<name>_markers.json.gz` |
| 3     | `differential_genes` | `--base_adata` (stage 2 h5ad) + `--diffgenes_contrasts` | per-contrast DE tables                            |

`--name` sets the identifier used in output file names; it defaults to the entry-point name.

> [!IMPORTANT]
> Running without `-entry` selects the upstream single-pass workflow, which the three-stage design
> replaced. It is retained for reference only and is no longer supported.

## Samplesheet input

The samplesheet is required by `-entry qc_clustering`. It is a comma-separated file with a header row
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
| `symbol_col`                    | Column holding gene symbols. Defaults to `index`. Two special values: `index` uses the matrix row names; `none` triggers gene symbol conversion via MyGene.info based on `geneid_col`. The values become the object index.                                                                              |
| `geneid_col`                    | Column holding gene identifiers. Defaults to `index`. Only used when `symbol_col` is `none`.                                                                                                                                                                                                            |
| `label_col`                     | Column holding cell type information. Defaults to `label`. If absent, the pipeline creates it and fills it with `unknown`. Renamed to `label` during execution.                                                                                                                                         |
| `unknown_label`                 | Value in `label_col` to treat as unknown. Defaults to `unknown`. Renamed to `unknown` during execution.                                                                                                                                                                                                 |
| `counts_layer`                  | Layer holding the raw counts matrix. Defaults to `X`.                                                                                                                                                                                                                                                   |
| `expected_cells`                | Expected number of cells, passed to CellBender for empty droplet detection.                                                                                                                                                                                                                             |
| `ambient_correction`            | Whether to run ambient RNA correction for this sample (`true` uses the globally configured method). Defaults to `true`.                                                                                                                                                                                 |
| `ambient_corrected_integration` | Whether to use ambient-corrected counts for integration for this sample, overriding the global `--ambient_corrected_integration`.                                                                                                                                                                       |
| `n_hvgs`                        | Number of highly variable genes for this sample. Defaults to `3000`.                                                                                                                                                                                                                                    |

> [!WARNING]
> **Per-sample QC filter columns are currently ignored.** The samplesheet accepts `min_genes`,
> `min_cells`, `min_counts_cell`, `min_counts_gene`, `max_mito_percentage` and
> `automatic_cell_filtering`, and they are read into the sample metadata, but the **global**
> `--min_*` / `--max_*` / `--automatic_cell_filtering` parameters are what actually reach the
> filtering step. Set your filters globally until this is fixed.

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
nextflow run UKDRI/scdownstream -r dev_ukdri -entry qc_clustering \
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
nextflow run UKDRI/scdownstream -r dev_ukdri -entry qc_clustering -profile apptainer -params-file params.yaml
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

On the UK DRI cluster, one `params_<entry>.yml` per stage is the standard pattern — see
[UK DRI usage](ukdri.md).

## `qc_clustering`

Stage 1: per-sample QC, cell type annotation, gene symbol unification, merging, scVI integration,
embeddings and Leiden clustering.

```bash
nextflow run UKDRI/scdownstream -r dev_ukdri -entry qc_clustering \
   -profile apptainer \
   --input samplesheet.csv \
   --name my_study \
   --species human \
   --outdir results/qc_clustering
```

Produces `results/qc_clustering/my_study_finalized.h5ad` (plus `.rds`), a Quarto QC/clustering report
under `report/`, and a MultiQC report.

| Parameter                         | Default          | Description                                                                                                                                                            |
| --------------------------------- | ---------------- | ---------------------------------------------------------------------------------------------------------------------------------------------------------------------- |
| `--input`                         | —                | Samplesheet (required).                                                                                                                                                |
| `--name`                          | `qc_clustering`  | Identifier used in output file names.                                                                                                                                  |
| `--species`                       | `human`          | `human` or `mouse`. Set explicitly for mouse data.                                                                                                                     |
| `--qc_only`                       | `false`          | Stop after per-sample QC and cell type annotation; skip merging, integration, embeddings, clustering and the Quarto report.                                            |
| `--mito_genes`                    | `null`           | File or pattern identifying mitochondrial genes, used for the mitochondrial percentage metric.                                                                         |
| `--ambient_correction`            | `decontx`        | `none`, `decontx`, `cellbender`, `soupx` or `scar`.                                                                                                                    |
| `--ambient_corrected_integration` | `false`          | Use ambient-corrected counts for integration rather than storing them as extra layers.                                                                                 |
| `--doublet_detection`             | `scrublet`       | Doublet detection method, or `none` to skip. See [Supported tool choices](#supported-tool-choices).                                                                    |
| `--automatic_cell_filtering`      | `false`          | Derive filtering thresholds automatically from N-MAD outlier detection instead of using the fixed thresholds below.                                                    |
| `--min_genes`                     | `200`            | Minimum genes per cell.                                                                                                                                                |
| `--max_genes`                     | `false`          | Maximum genes per cell (`false` disables).                                                                                                                             |
| `--min_cells`                     | `5`              | Minimum cells per gene.                                                                                                                                                |
| `--min_counts`                    | `0`              | Minimum counts per cell.                                                                                                                                               |
| `--max_counts`                    | `false`          | Maximum counts per cell (`false` disables).                                                                                                                            |
| `--min_counts_gene`               | `0`              | Minimum counts per gene.                                                                                                                                               |
| `--max_mito_percentage`           | `25`             | Maximum percentage of mitochondrial counts per cell.                                                                                                                   |
| `--symbol_col`                    | `index`          | Default gene symbol column for samples that do not set their own.                                                                                                      |
| `--unify_gene_symbols`            | `false`          | Unify gene symbols across samples using HUGO.                                                                                                                          |
| `--duplicate_var_resolution`      | `sum`            | How to resolve duplicate gene symbols: `mean`, `sum`, `max` or `make_unique`.                                                                                          |
| `--aggregate_isoforms`            | `false`          | Aggregate isoform-level features.                                                                                                                                      |
| `--n_hvgs`                        | `3000`           | Highly variable genes used for the PCA/UMAP embedding.                                                                                                                 |
| `--integration_methods`           | `scvi`           | Integration method. Keep `scvi` — see [Supported tool choices](#supported-tool-choices).                                                                               |
| `--integration_hvgs`              | `5000`           | Highly variable genes used for integration.                                                                                                                            |
| `--clustering_resolutions`        | `0.5,1.0`        | Comma-separated Leiden resolutions. One `leiden_<res>` column is written per resolution, and the **first** resolution is copied to `leiden` as the default clustering. |
| `--cluster_neighbors`             | `neighbors_scvi` | Neighbour graph used for clustering. The default assumes scVI integration.                                                                                             |
| `--celltypist_model`              | `null`           | Comma-separated [CellTypist](https://www.celltypist.org/models) model names, or a path to a custom `.pkl` model.                                                       |
| `--celldex_reference`             | `null`           | CSV describing celldex references for SingleR — see [Cell type annotation](#cell-type-annotation).                                                                     |
| `--scvi_model` / `--scanvi_model` | `null`           | Pre-trained model for reference mapping — see [Reference mapping](#reference-mapping).                                                                                 |
| `--save_intermediates`            | `false`          | Publish intermediate objects as well as the final ones.                                                                                                                |

scVI itself is tuned with `--scvi_n_latent` (30), `--scvi_n_hidden` (128), `--scvi_n_layers` (2),
`--scvi_dispersion` (`gene`), `--scvi_gene_likelihood` (`zinb`), `--scvi_max_epochs`,
`--scvi_categorical_covariates` and `--scvi_continuous_covariates`.

> [!NOTE]
> `--species` defaults to `human`. Set it explicitly for mouse data: the human default is applied
> silently, and gene set enrichment and LIANA+ ortholog mapping would then use the wrong species.

## `downstream`

Stage 2: marker genes per cluster, gene set enrichment, and LIANA+ cell–cell communication.

```bash
nextflow run UKDRI/scdownstream -r dev_ukdri -entry downstream \
   -profile apptainer \
   --base_adata results/qc_clustering/my_study_finalized.h5ad \
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
| `--ortholog_hcop_directory`       | `/nfsdata/genome/hcop/` | Directory of HCOP ortholog tables for LIANA+ — see [Reference data](#reference-data).                                                                   |
| `--markers_uns_key`               | `rank_genes_groups`     | `uns` key holding the marker results to export.                                                                                                         |
| `--markers_thr_adj_pvalue`        | `0.05`                  | Adjusted p-value threshold for the exported markers.                                                                                                    |
| `--markers_n_top`                 | `100`                   | Number of top markers per group to export.                                                                                                              |
| `--markers_pct_nz`                | `0.1`                   | Minimum fraction of non-zero expressing cells for an exported marker.                                                                                   |
| `--markers_min_logfc`             | `0`                     | Minimum log fold change for an exported marker.                                                                                                         |
| `--report_table_row_limit`        | `250`                   | Maximum rows shown per table in the HTML report.                                                                                                        |

Marker genes are computed with `scanpy.tl.rank_genes_groups` using the **Wilcoxon** test.

## `differential_genes`

Stage 3: pseudobulk aggregation with decoupler, then PyDESeq2 differential expression for every
combination of group label and contrast.

```bash
nextflow run UKDRI/scdownstream -r dev_ukdri -entry differential_genes \
   -profile apptainer \
   --base_adata results/downstream/my_study_finalized.h5ad \
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
| `exclude_samples_col`    | no       | Accepted but **currently inactive**.                                                                                            |
| `exclude_samples_values` | no       | Accepted but **currently inactive**.                                                                                            |

The design formula is `~ variable [+ blocking…]`, fitted with `refit_cooks=True`.

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
| Ambient RNA correction | `--ambient_correction`                      | `decontx` (default), `soupx`, `cellbender`, `scar`, `none` |
| Integration            | `--integration_methods`                     | `scvi` (default)                                           |
| Clustering             | `--clustering_resolutions`                  | Leiden, one or more resolutions                            |
| Cell type annotation   | `--celltypist_model`, `--celldex_reference` | CellTypist, SingleR / celldex                              |

Two consequences of the current tool set are worth knowing before you interpret results:

- **Doublets are annotated, not removed.** scrublet writes its scores and its `predicted_doublet`
  call into the object; no cells are dropped. Filter on that annotation yourself in downstream
  analysis. `--doublet_detection_threshold` has no effect.
- **scVI drives everything after the merge.** The embeddings, clustering and reports are built on the
  scVI latent space, and `--cluster_neighbors` defaults to `neighbors_scvi` accordingly. Keep `scvi`
  in `--integration_methods`; if it is omitted, the post-merge steps produce no output and no error.

**Inactive parameters.** The following are accepted but currently have no effect:
`--doublet_detection_threshold`, `--skip_enrichment`, `--skip_liana`, `--skip_rankgenesgroups`,
`--pseudobulk`, `--pseudobulk_groupby_labels`, `--pseudobulk_min_num_cells`, `--cluster_per_label`,
`--cluster_global`, `--base_embeddings`, `--base_label_col`. `--prep_cellxgene` is no longer
supported and should be left at its default.

## Local container images

Several modules added by this fork use custom container images that are not on a public registry.
The `--singularity_cache_dir` parameter tells them where to find locally built `.sif` files:

```bash
--singularity_cache_dir /nfsdata/apptainer
```

It defaults to `$NXF_SINGULARITY_CACHEDIR`, or `$NXF_APPTAINER_CACHEDIR` if that is unset, so
exporting either environment variable is usually enough. When set and running under the
`singularity` or `apptainer` profile, modules use `<singularity_cache_dir>/<image>.sif` instead of
pulling a remote container.

The images referenced are:

| Image                           | Used by                                     |
| ------------------------------- | ------------------------------------------- |
| `scanpy_1.11.4_coreinf_0.3.sif` | the Quarto report modules                   |
| `scanpy_1.11.4_coreinf_0.1.sif` | `SCANPY_ENRICH`, `SCANPY_EXPORT_MARKERS`    |
| `decoupler_latest.sif`          | `DECOUPLER_PSEUDOBULK`, `FILTER_PSEUDOBULK` |
| `pydeseq2_latest.sif`           | `DIFFERENTIAL_GENES_PER_CONTRAST`           |

> [!WARNING]
> The decoupler and PyDESeq2 modules fall back to a public registry image if no local `.sif` is
> found, but **five processes do not**: `SCANPY_GENERATE_REPORT`, `SCANPY_GENERATE_REPORT_QC`,
> `PYDESEQ2_GENERATE_REPORT`, `SCANPY_ENRICH` and `SCANPY_EXPORT_MARKERS` fall back to hard-coded
> absolute paths on the UK DRI filesystem and ignore the container engine in use. Off the UK DRI
> cluster you must build these images and point `--singularity_cache_dir` at a directory holding them
> under exactly the names above.

## Reference data

**HCOP orthologs (LIANA+).** `--ortholog_hcop_directory` defaults to the UK DRI path
`/nfsdata/genome/hcop/`. LIANA+ reads `<directory>/human_<species>_hcop_fifteen_column.txt.gz` from
it to map its human-derived ligand–receptor resource onto non-human data. Off-site runs must
override the parameter and provide the corresponding
[HCOP](https://www.genenames.org/tools/hcop/) table.

**CellTypist models** are downloaded at runtime unless `--celltypist_model` is given a local `.pkl`
path.

**celldex references** may be given by name (downloaded at runtime) or as paths to pre-downloaded tar
archives — useful where the compute nodes have no internet access.

## Cell type annotation

### CellTypist

Specify one or more models with `--celltypist_model` (comma-separated). Available models are listed
at [celltypist.org](https://www.celltypist.org/models). A path containing `/` is treated as a custom
model file.

### SingleR

`--celldex_reference` takes a CSV describing the celldex references to use. The available references
are described in the
[celldex package documentation](https://bioconductor.org/packages/devel/data/experiment/manuals/celldex/man/celldex.pdf).

Referring to references by name:

```csv title="celldex_references.csv"
id,label,reference,version
hpca,label.main,hpca,2024-02-26
monaco_immune,label.fine,monaco_immune,2024-02-26
```

Referring to pre-downloaded tar archives by path:

```csv title="celldex_references.csv"
hpca,label.main,/path/to/hpca.tar
monaco_immune,label.fine,/path/to/monaco_immune.tar
```

> [!NOTE]
> Cell type predictions are merged into the per-sample objects as `obs` columns during the
> finalisation step of `qc_clustering`.

## Ambient RNA correction

Ambient RNA correction removes contaminating RNA from cell-free droplets. Select the method globally
with `--ambient_correction` (`decontx` by default; also `cellbender`, `soupx`, `scar`, or `none`):

```bash
nextflow run UKDRI/scdownstream -r dev_ukdri -entry qc_clustering \
   --ambient_correction cellbender --input samplesheet.csv --outdir results
```

Correction can be disabled per sample from the samplesheet:

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

Selects the stage to run: `qc_clustering`, `downstream` or `differential_genes`. Always required —
see [Choosing an entry point](#choosing-an-entry-point).

### `-profile`

Chooses a configuration profile. Profiles bundled with the pipeline select how software is provided:

- `apptainer` — [Apptainer](https://apptainer.org/) (the UK DRI default)
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

For the custom images this fork relies on, prefer `--singularity_cache_dir` — see
[Local container images](#local-container-images).

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
> such as `--singularity_cache_dir` and `--ortholog_hcop_directory`.

## Running in the background

Nextflow must keep running until the pipeline finishes. Use `-bg` to detach it from your terminal
(logs go to a file), or run it inside `screen`/`tmux`. On an HPC, run Nextflow itself as a cluster job
from which it submits the individual tasks — this is the UK DRI pattern, see [UK DRI usage](ukdri.md).

## Nextflow memory requirements

The Nextflow JVM can request a large amount of memory. Limiting it in your shell profile (typically
`~/.bashrc`) is recommended:

```bash
NXF_OPTS='-Xms1g -Xmx4g'
```
