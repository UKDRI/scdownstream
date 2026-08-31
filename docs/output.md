# UK DRI scdownstream: Output

## Introduction

This document describes the output produced by the pipeline. All paths are relative to the directory
given by `--outdir`.

Because the pipeline is split into [three entry points](usage.md#choosing-an-entry-point), each stage
writes into its own results directory and produces a different subset of the tree below. The sections
are grouped by stage.

> [!IMPORTANT]
> Most intermediate objects are only published when `--save_intermediates` is set. Without it you get
> the finalized objects, the HTML reports, MultiQC and `pipeline_info/` — which is usually what you
> want.

## Stage 1 — `qc_clustering`

Produces, at the top level of the results directory:

<details markdown="1">
<summary>Output files</summary>

- `<name>_finalized.h5ad`: the integrated, clustered AnnData object. **This is the input to
  `-entry downstream`.**
- `<name>_finalized.rds`: a SingleCellExperiment version of the same object.

</details>

`<name>` is the value of `--name` (defaulting to `qc_clustering`).

> [!NOTE]
> RDS conversion runs with `errorStrategy 'ignore'`, so the `.rds` file is best-effort: if conversion
> fails the run continues and only the `.h5ad` is produced.

### Input loading

<details markdown="1">
<summary>Output files</summary>

- `load_h5ad/`: the result of converting RDS, 10x h5 and CSV inputs to h5ad.

</details>

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
  - `custom_thresholds/`: the result of cell and gene filtering, plus QC plots. When
    `--automatic_cell_filtering` is set, the N-MAD-derived thresholds used are reported here.
  - `qc_filtered/`: QC metrics and plots after filtering, also fed into MultiQC.
  - `finalized/`: the per-sample objects with cell type predictions merged back in.

</details>

Note the step order in this fork: **doublet detection runs before ambient RNA correction and
filtering**.

> [!WARNING]
> Doublets are **annotated, not removed**. The `obs` columns written by scrublet mark predicted
> doublets, but no cells are dropped — see
> [Supported tool choices](usage.md#supported-tool-choices).

### Cell type annotation

<details markdown="1">
<summary>Output files</summary>

- `celltypes/`
  - `celltypist/`: the annotated `h5ad` and the predictions as a `pkl`.
  - `singler/`: the annotated `h5ad`, the predictions as a `csv`, and diagnostic plots as `pdf`.
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

### Reports

<details markdown="1">
<summary>Output files</summary>

- `report/<name>_scdownstream_report.html`: a self-contained [Quarto](https://quarto.org/) report
  covering per-sample QC, filtering, integration, the UMAP embedding and the clustering at each
  resolution. Tables are searchable and capped at `--report_table_row_limit` rows.

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

- `<name>_finalized.h5ad` / `<name>_finalized.rds`: the object with marker genes, enrichment results
  and LIANA+ output added. **This is the input to `-entry differential_genes`.**
- `<name>_markers.json.gz`: the marker genes per group, filtered by `--markers_thr_adj_pvalue`,
  `--markers_n_top`, `--markers_pct_nz` and `--markers_min_logfc`, in a compact JSON form for
  downstream tooling.
- `scanpy/<name>/`: the object after `rank_genes_groups` (Wilcoxon) and after gene set enrichment.
- `per_group/<name>/liana/`: [LIANA+](https://liana-py.readthedocs.io/) rank-aggregate cell–cell
  communication results.
- `report/<name>_scdownstream_report.html`: a Quarto report covering the marker genes per group, the
  enrichment results and the cell–cell communication analysis.

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
  every group label — pseudobulk composition, per-contrast result tables and summary counts.

</details>

Only the `pseudobulk/` intermediates are gated on `--save_intermediates`; the DE tables and the
report are always published.

> [!WARNING]
> A (group label, contrast) combination with fewer than `--diffgenes_min_samples` pseudobulk samples
> on either side is **skipped silently** and produces no `.tsv`. The report is driven by a
> `de_manifest.tsv` mapping each result file to its group label and contrast, so the set of tables in
> the report is the authoritative record of what was actually tested. Compare it against the contrasts
> you supplied.

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
The three entry points do not write them: `cluster_dimred/`, `pseudobulking/`,
`per_group/<id>/paga/`, `finalized/`, and the
`quality_control/doublet_detection/{solo,doubletdetection,scds}/` directories.
