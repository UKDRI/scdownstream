# UK DRI scdownstream: Output

## Introduction

This document describes the output produced by the pipeline. All paths are relative to the directory
given by `--outdir`.

Because the pipeline is split into [three stages](usage.md#choosing-a-stage), each stage
writes into its own results directory and produces a different subset of the tree below. The sections
are grouped by stage.

> [!IMPORTANT]
> Most intermediate objects are only published when `--save_intermediates` is set. Without it you get
> the finalized objects, the HTML reports, MultiQC and `pipeline_info/` — which is usually what you
> want. Intermediate h5ad files are written uncompressed and may hold dense matrices; only the
> finalized `<name>_qc_clustering.h5ad` and `<name>_downstream.h5ad` are guaranteed sparse and
> gzip-compressed.

## Stage 1 — `qc_clustering`

Produces, at the top level of the results directory:

<details markdown="1">
<summary>Output files</summary>

- `<name>_qc_clustering.h5ad`: the integrated, clustered AnnData object, gzip-compressed, with `X`,
  every layer and `raw.X` stored as sparse CSR matrices. **This is the input to `--stage downstream`.**
- `<name>_qc_clustering.rds`: a SingleCellExperiment version of the same object.

</details>

`<name>` is the value of `--name` (`scdownstream` if it is not set); the same applies to stage 2.

> [!NOTE]
> RDS conversion runs with `errorStrategy 'ignore'`, so the `.rds` file is best-effort: if conversion
> fails the run continues and only the `.h5ad` is produced.

### Input loading

<details markdown="1">
<summary>Output files</summary>

- `load_h5ad/`: the result of converting RDS, 10x h5 and CSV inputs to h5ad.
- `load/`: the per-sample objects with the `sample` obs column added (`<sample>_sample.h5ad`) and,
  with `--metadata`, the per-sample metadata columns added (`<sample>_metadata.h5ad`).

</details>

The `sample` column and the `--metadata` columns are added to `obs` before QC. So they are in every
later object, including `<name>_qc_clustering.h5ad` and `<name>_downstream.h5ad`. See
[Per-sample metadata](usage.md#per-sample-metadata).

### Quality control

<details markdown="1">
<summary>Output files</summary>

- `quality_control/`
  - `sizes/`: collected per-sample object sizes.
  - `empty_droplet_removal/`: CellBender empty droplet detection. Only produced when no `filtered`
    matrix was given in the samplesheet.
  - `qc_raw/`: QC metrics and plots for the raw input data, also fed into MultiQC.
  - `doublet_detection/scrublet/`: scrublet output — the annotated `h5ad`, the per-cell doublet
    scores and predictions, and a `*_mqc.json` carrying the score distribution plot for MultiQC.
  - `ambient_rna_removal/`: ambient RNA correction results (decontX, SoupX, CellBender or scAR
    depending on `--ambient_correction`).
  - `custom_thresholds/`: the result of cell and gene filtering, plus QC plots. The thresholds
    applied are reported here and in MultiQC; with `--automatic_cell_filtering` these are the N-MAD
    thresholds after the `--min_genes` / `--min_counts` floors, and the report says when a floor
    was used. The number of cells removed because they had no counts left after gene filtering is
    reported too.
  - `qc_filtered/`: QC metrics and plots after filtering, also fed into MultiQC.
  - `finalized/`: the per-sample objects with cell type predictions merged back in.

</details>

Note the step order in this fork: **doublet detection runs before ambient RNA correction and
filtering**.

> [!WARNING]
> By default, doublets are **annotated, not removed**: the `obs` columns written by scrublet
> (`doublet_score`, `predicted_doublet`) mark predicted doublets, but no cells are dropped. With
> `--doublet_removal` they are removed right after detection. With `--filtering_keep_outliers`,
> cells that fail the QC thresholds are kept and marked in a bool `obs["outlier"]` column instead
> of removed. See [Doublets and QC outliers](usage.md#doublets-and-qc-outliers).

### Cell type annotation

<details markdown="1">
<summary>Output files</summary>

- `celltypes/`
  - `celltypist/`: the annotated `h5ad` and the predictions as a `pkl`.
- `adata/`: objects produced by merging predictions back into the sample AnnData.

</details>

### Column harmonisation

<details markdown="1">
<summary>Output files</summary>

- `unify/`
  - `*.h5ad`: the per-sample objects after their gene symbol, batch and cell type columns have been
    harmonised and duplicate gene symbols resolved, ready for merging.

</details>

Gene symbols are resolved from each sample's `symbol_col` / `geneid_col` (via MyGene.info where
`symbol_col` is `none`), duplicates are combined according to `--duplicate_var_resolution`, and
isoforms are optionally aggregated with `--aggregate_isoforms`.

### Merging and integration

<details markdown="1">
<summary>Output files</summary>

- `combine/`
  - `merge/`
    - `merged_inner.h5ad`: samples merged on the intersection of genes. Used for integration.
    - `merged_outer.h5ad`: samples merged over all genes. Used as the base for the final object.
    - `merged_sample_genes.png`: UpSet plot of gene overlap between samples.
  - `merge_emb/<id>/`: merged embeddings, when embeddings from multiple sources are combined.
  - `integrate/`
    - `input_hvg/`: the highly-variable-gene subset passed to the integration tools, as `h5ad` and
      `rds`.
    - `scvi/`: the scVI-integrated object and `X_scvi.pkl`, the latent representation.

</details>

scVI is the curated integration method — see
[Supported tool choices](usage.md#supported-tool-choices).

### Embeddings and clustering

<details markdown="1">
<summary>Output files</summary>

- `scanpy/<name>/`: the object after each embedding and clustering step — log-normalisation, HVG
  selection, PCA (including loadings), the neighbour graph, UMAP, and Leiden clustering.

</details>

Leiden clustering writes one `leiden_<resolution>` column per entry in
`--clustering_resolutions`, and copies the **first** resolution into `leiden` as the default
clustering. `--selected_clustering` in stage 2 picks which of these to use.

Embeddings stored in the object's `obsm`:

| Key           | Content                                                                                   |
| ------------- | ----------------------------------------------------------------------------------------- |
| `X_scvi`      | scVI latent representation (also as `X_emb`).                                             |
| `X_umap_scvi` | UMAP on the scVI neighbour graph, the graph Leiden clustering uses. Plotted by default.   |
| `X_pca`       | PCA embedding.                                                                            |
| `X_umap_pca`  | UMAP on the PCA neighbour graph.                                                          |
| `X_umap`      | Copy of `X_umap_pca`, kept only because tools such as `sc.pl.umap` expect this key.      |

The PCA-based UMAP is stored as `X_umap_pca` so its name says what it is built on; `X_umap` alone
does not, and without checking this pipeline's steps you would not know what it was computed from.

### Reports

<details markdown="1">
<summary>Output files</summary>

- `report/<name>_qc_clustering_report.html`: a self-contained [Quarto](https://quarto.org/) report
  covering per-sample QC, filtering, integration, the UMAP embedding and the clustering at each
  resolution. The PCA section has an elbow plot of the variance explained per PC, marking the
  number of PCs the PCA graph used (`--neighbors_n_pcs`). UMAP plots use `--umap_for_plots`
  (default `X_umap_scvi`). With `--filtering_keep_outliers`, a "QC outliers" part of the UMAP
  section shows where the marked outliers fall, how many there are per sample, and their QC
  metrics compared with the other cells. With `--umap_color_by`,
  a "Requested UMAP plots" section plots those `obs` columns on each UMAP in
  `--umap_color_by_embeddings` (by default before and after integration). Tables are searchable
  and capped at `--report_table_row_limit` rows.

</details>

### MultiQC

<details markdown="1">
<summary>Output files</summary>

- `multiqc/`
  - `multiqc_report.html`: a standalone HTML report viewable in a browser.
  - `multiqc_data/`: parsed statistics from the tools used.
  - `multiqc_plots/`: static images from the report.

</details>

[MultiQC](http://multiqc.info) collates per-sample QC across the pipeline — raw and filtered cell and
gene counts, mitochondrial content, doublet score distributions — into one report, and also records
the software versions used. Coverage is partial: the Quarto report above is the more complete view of
a `qc_clustering` run.

## Stage 2 — `downstream`

<details markdown="1">
<summary>Output files</summary>

- `<name>_downstream.h5ad` (gzip-compressed, sparse CSR `X`, layers and `raw.X`) / `<name>_downstream.rds`: the
  object with marker genes, enrichment results and LIANA+ output added. **This is the input to `--stage differential_genes`.**
- `<name>_downstream_markers.json.gz`: the marker genes per group (one set per cluster, id
  `genes--cluster--<cluster>`; the separator is `--markers_sep`), filtered by `--markers_thr_adj_pvalue`,
  `--markers_n_top`, `--markers_pct_nz` and `--markers_min_logfc`, in a compact JSON form for
  downstream tooling.
- `scanpy/<name>/`: the object after `rank_genes_groups` (Wilcoxon) and after gene set enrichment.
- `per_group/<name>/liana/`: [LIANA+](https://liana-py.readthedocs.io/) rank-aggregate cell–cell
  communication results.
- `report/<name>_downstream_report.html`: a Quarto report covering the marker genes per group, the
  enrichment results and the cell–cell communication analysis. With `--umap_color_by`, it also
  has the "Requested UMAP plots" section described for stage 1.

</details>

Marker genes and enrichment are computed over the grouping named by `--selected_clustering`, so a
run grouped by cell type and a run grouped by `leiden_1.0` should be written to different `--outdir`
directories (or given different `--name` values) to avoid overwriting each other.

## Stage 3 — `differential_genes`

<details markdown="1">
<summary>Output files</summary>

- `differential_genes/`
  - `pseudobulk/`
    - `<name>_pseudobulk.h5ad`: raw pseudobulk profiles, one observation per sample × group label,
      summed from the `counts` layer.
    - `<name>_pseudobulk_filtered.h5ad`: the same after dropping profiles below
      `--diffgenes_min_counts` / `--diffgenes_min_cells`.
    - `per_group/`: one h5ad per distinct value of `--diffgenes_group_col`.
  - `deseq2/<group_label>/`
    - `<name>_<group_label>_<contrast>.tsv`: the PyDESeq2 result table for that group label and
      contrast — one row per gene, with columns `feature`, `baseMean`, `log2FoldChange`, `lfcSE`,
      `stat`, `pvalue`, `padj`. Log fold changes are **target vs reference**.
- `report/<name>_differential_genes_report.html`: a Quarto report summarising every contrast across
  every group label. In order:
  - **Contrasts**: one row per contrast with the variable, target and reference, the direction of
    the log2 fold change in words (e.g. "> 0: higher in AD than in ctrl"), blocking variables and
    excluded samples.
  - **Significant genes overview**: group label × contrast, with the number of significant up- and
    down-regulated genes, or "skipped" where the contrast was not tested.
  - **Pseudobulk summary**: the profiles per sample, group label and contrast variable.
  - **Volcano plots**: one section per group label, headed `<group column> = <value>` (e.g.
    `celltypist:majority_voting = Microglia`). Each starts with a table of significant up / down /
    tested genes per contrast, including skipped contrasts, followed by each contrast's volcano plot
    and its up- and down-regulated gene tables. A volcano plot with no significant genes is
    collapsed; click to show it.

  "Significant" means `padj` < 0.05 and |log2FoldChange| > 1, the same rule for the counts, the
  volcano plots and the gene tables.
- `<name>_differential_genes.json.gz`: the top differential genes per group label and
  contrast, in the same compact JSON form as `<name>_downstream_markers.json.gz`. Each group label ×
  contrast gives two sets: `genes--<contrast>--<group_label>--up` (higher in the target group,
  ranked by the Wald statistic) and `genes--<contrast>--<group_label>--down` (lower in the target
  group), e.g. `genes--AD_vs_ctrl--Microglia--up`. The separator (`--` by default, because contrast ids often
  contain `__`) is
  `--markers_sep` and is recorded in `params.set_id_separator`, so the ids can be split into their
  parts; the run stops if a group value or contrast id contains it. Genes need an adjusted p-value
  below `--diffgenes_markers_thr_adj_pvalue` and an absolute log2 fold change of at least
  `--diffgenes_markers_min_logfc`; each set holds at most `--diffgenes_markers_n_top` genes, and
  empty sets are left out. Per gene: `list_name`, `list_p_value`, `list_p_value_adjusted`,
  `list_logFC` (log2 fold change), `list_score` (Wald statistic) and `list_base_mean` (pseudobulk
  has no per-cell detection rates, so there are no `pct_nz` lists).

</details>

Only the `pseudobulk/` intermediates are gated on `--save_intermediates`; the DE tables, the
marker JSON and the report are always published.

> [!WARNING]
> A (group label, contrast) combination with fewer than `--diffgenes_min_samples` pseudobulk samples
> on either side is **skipped** and produces no `.tsv`; the run does not fail. The report lists every
> contrast of the contrasts file for every group label and marks such combinations as "skipped" in the
> overview and in the group's table, so you can see what was not tested.

## Pipeline information

<details markdown="1">
<summary>Output files</summary>

- `pipeline_info/`
  - Reports generated by Nextflow, each suffixed with the run timestamp:
    `execution_report_*.html`, `execution_timeline_*.html`, `execution_trace_*.txt` and
    `pipeline_dag_*.html`.
  - Software versions: `nf_core_scdownstream_software_mqc_versions.yml` (`qc_clustering`) or
    `nf_core_scdownstream_differential_genes_software_mqc_versions.yml` (`differential_genes`).
  - `pipeline_report.html` / `pipeline_report.txt`: only produced by `qc_clustering`, and only when
    `--email` / `--email_on_fail` is used.
  - The validated samplesheet: `samplesheet.valid.csv`.
  - The parameters used for the run: `params.json`.

</details>

[Nextflow](https://www.nextflow.io/docs/latest/tracing.html) generates these reports to help
troubleshoot failures and to record launch commands, run times and resource usage.

## Directories referenced elsewhere but not produced

Some publishing rules inherited from upstream target steps that are not part of the curated workflow.
The three stages do not write them: `cluster_dimred/`, `pseudobulking/`,
`per_group/<id>/paga/`, `finalized/`, and the
`quality_control/doublet_detection/{solo,doubletdetection,scds}/` directories.
