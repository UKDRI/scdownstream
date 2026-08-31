# UK DRI scdownstream

_A UK DRI fork of [nf-core/scdownstream](https://github.com/nf-core/scdownstream)._

> [!IMPORTANT]
> This pipeline began life as a fork of [nf-core/scdownstream](https://github.com/nf-core/scdownstream) and has
> since diverged substantially. It is now developed independently by **UK DRI Informatics** and follows its own release and configuration conventions.
>
> Please raise issues and questions in [this repository](https://github.com/UKDRI/scdownstream). For
> the original, nf-core community-supported pipeline, see
> [nf-co.re/scdownstream](https://nf-co.re/scdownstream).
>
> The original authors and contributors are credited in
> [`ACKNOWLEDGEMENTS.md`](ACKNOWLEDGEMENTS.md).

[![Nextflow](https://img.shields.io/badge/nextflow-%E2%89%A524.10.5-green?style=flat&logo=nextflow&logoColor=white&color=%230DC09D)](https://www.nextflow.io/)
[![run with apptainer](https://img.shields.io/badge/run%20with-apptainer-1d355c.svg?labelColor=000000)](https://apptainer.org/)
[![License: MIT](https://img.shields.io/badge/license-MIT-blue.svg?labelColor=000000)](LICENSE)

## Introduction

**UK DRI scdownstream** is a Nextflow pipeline for the downstream analysis of processed single-cell
RNA-seq data. It takes per-sample count matrices (h5ad, 10x h5, SingleCellExperiment / Seurat RDS, or
CSV) and carries them through quality control, ambient RNA correction, doublet detection, cell type
annotation, scVI integration, dimensionality reduction and clustering, then on through marker genes,
gene set enrichment and cell–cell communication, and finally pseudobulk differential expression.

The input is ideally the output of [nf-core/scrnaseq](https://nf-co.re/scrnaseq) or a similar
per-sample count matrix object files. We recommend supplying it as AnnData (`.h5ad`) objects: AnnData is
part of the [scverse](https://scverse.org/) ecosystem that this pipeline is built on, and it scales
to large count matrices without the size limits that SingleCellExperiment and Seurat objects run into
in R.

The pipeline is organised as **three sequential stages**, each its own Nextflow entry point, so that
long-running analyses can be checkpointed and re-run independently.

It inherits the learnings and implementations of the following pipelines (alphabetical), via
nf-core/scdownstream:

- [panpipes](https://github.com/DendrouLab/panpipes)
- [scFlow](https://combiz.github.io/scFlow/)
- [scRAFIKI](https://github.com/Mye-InfoBank/scRAFIKI)
- [YASCP](https://github.com/wtsi-hgi/yascp)

## Entry points

The pipeline is always invoked with an explicit `-entry`. The three stages are **sequential**: each
consumes the previous stage's finalized `.h5ad` via `--base_adata`.

| Stage | Entry point                 | Required input                                                        | Main output                                                             |
| ----- | --------------------------- | --------------------------------------------------------------------- | ----------------------------------------------------------------------- |
| 1     | `-entry qc_clustering`      | `--input samplesheet.csv`                                             | `<name>_finalized.h5ad`, MultiQC report, QC/clustering HTML report      |
| 2     | `-entry downstream`         | `--base_adata <stage 1 h5ad>`                                         | `<name>_finalized.h5ad`, `<name>_markers.json.gz`, analysis HTML report |
| 3     | `-entry differential_genes` | `--base_adata <stage 2 h5ad>` + `--diffgenes_contrasts contrasts.tsv` | per-contrast DE tables, pseudobulk h5ad, DE HTML report                 |

```text
samplesheet.csv
      │
      ▼  -entry qc_clustering
<name>_finalized.h5ad  ── QC, integration (scVI), UMAP, Leiden clustering
      │
      ▼  -entry downstream        --base_adata
<name>_finalized.h5ad  ── marker genes, enrichment, LIANA+ cell–cell communication
      │
      ▼  -entry differential_genes  --base_adata  --diffgenes_contrasts
differential_genes/    ── decoupler pseudobulk → PyDESeq2 per group × contrast
```

> [!IMPORTANT]
> **Always pass `-entry`.** Running `nextflow run …` without it selects the upstream single-pass
> workflow ([`workflows/scdownstream.nf`](workflows/scdownstream.nf)), which the three-stage design
> replaced. It is retained for reference only and is no longer supported.

## Pipeline steps

### Stage 1 — `qc_clustering`

1. Load and convert inputs to h5ad (h5ad, 10x h5, RDS, CSV)
2. Per-sample quality control
   1. QC metrics for raw counts ([`MultiQC`](http://multiqc.info/))
   2. Doublet detection — [scrublet](https://scanpy.readthedocs.io/en/stable/api/generated/scanpy.pp.scrublet.html)
      (doublets are **annotated, not removed** — see [Status](#status-and-known-limitations))
   3. Ambient RNA correction — [decontX](https://bioconductor.org/packages/release/bioc/html/decontX.html)
      (default), [soupX](https://cran.r-project.org/web/packages/SoupX/readme/README.html),
      [CellBender](https://cellbender.readthedocs.io/en/latest/),
      [scAR](https://docs.scvi-tools.org/en/stable/user_guide/models/scar.html)
   4. Cell filtering — fixed thresholds, or automatic N-MAD outlier thresholds via
      `--automatic_cell_filtering`
3. Cell type annotation — [CellTypist](https://www.celltypist.org/) and/or
   [SingleR](https://bioconductor.org/packages/release/bioc/html/SingleR.html) with
   [celldex](https://bioconductor.org/packages/release/data/experiment/html/celldex.html) references
4. Harmonise gene symbol, batch and cell type columns across samples, resolving duplicate gene
   symbols
5. Merge and integrate — [scVI](https://docs.scvi-tools.org/en/stable/user_guide/models/scvi.html)
6. Embeddings and clustering — log-normalisation, HVG selection, PCA, neighbours,
   [UMAP](https://scanpy.readthedocs.io/en/stable/generated/scanpy.tl.umap.html),
   [Leiden](https://scanpy.readthedocs.io/en/stable/generated/scanpy.tl.leiden.html) at every
   resolution in `--clustering_resolutions`
7. Reports — a [Quarto](https://quarto.org/) QC/clustering report and a
   [`MultiQC`](http://multiqc.info/) report

`--qc_only` stops after step 3, producing per-sample QC objects without merging or integration.

### Stage 2 — `downstream`

1. Marker genes per cluster — `scanpy.tl.rank_genes_groups` (Wilcoxon) on the clustering named by
   `--selected_clustering`
2. Gene set enrichment over those markers
3. Cell–cell communication — [LIANA+](https://liana-py.readthedocs.io/) rank aggregation, using HCOP
   orthologs for non-human species
4. Marker export to `<name>_markers.json.gz`, and h5ad/RDS finalisation
5. A Quarto analysis report

### Stage 3 — `differential_genes`

1. Pseudobulk aggregation per sample × group label — [decoupler](https://decoupler-py.readthedocs.io/)
   (`counts` layer, sum mode)
2. Drop pseudobulk samples below `--diffgenes_min_counts` / `--diffgenes_min_cells`
3. Split into one object per group label (cluster or cell type, `--diffgenes_group_col`)
4. Differential expression per group label × contrast —
   [PyDESeq2](https://pydeseq2.readthedocs.io/), design `~ variable [+ blocking…]`
5. A Quarto differential expression report across all results

Contrasts are supplied as a TSV following the
[nf-core/differentialabundance](https://nf-co.re/differentialabundance) definition — see
[`assets/contrasts.tsv`](assets/contrasts.tsv) and the
[usage documentation](docs/usage.md#differential_genes).

## Tool choices

This fork focuses on a curated set of tools — the approaches we have validated on UK DRI data:

| Step                    | Supported                                                     |
| ----------------------- | ------------------------------------------------------------- |
| Doublet detection       | `scrublet`                                                    |
| Ambient RNA correction  | `decontx` (default), `soupx`, `cellbender`, `scar`, or `none` |
| Integration             | `scvi`                                                        |
| Clustering              | Leiden, at multiple resolutions                               |
| Cell type annotation    | CellTypist, SingleR / celldex                                 |
| Differential expression | PyDESeq2 on decoupler pseudobulk                              |

Further tools will be added as they are curated and validated. Until then, please use the values
above — see [Supported tool choices](docs/usage.md#supported-tool-choices) for the details.

## Quick start

> [!NOTE]
> If you are new to Nextflow, see the [Nextflow documentation](https://www.nextflow.io/docs/latest/)
> for installation. If you are unsure about the `filtered` / `unfiltered` distinction in the
> samplesheet, see [Filtered and unfiltered matrices](docs/usage.md#filtered-and-unfiltered-matrices).

Prepare a samplesheet describing your per-sample matrices:

```csv title="samplesheet.csv"
sample,unfiltered
sample1,/absolute/path/to/sample1.h5ad
sample2,/absolute/path/to/sample2.h5
sample3,relative/path/to/sample3.rds
sample4,/absolute/path/to/sample4.csv
```

Each entry is an h5ad, 10x h5, RDS or CSV file. RDS files may contain any object convertible to a
SingleCellExperiment via [Seurat's `as.SingleCellExperiment`](https://satijalab.org/seurat/reference/as.singlecellexperiment).
CSV files should hold a matrix with genes as columns and cells as rows, the first column being cell
barcodes.

**Stage 1 — QC, integration and clustering:**

```bash
nextflow run UKDRI/scdownstream -r dev_ukdri -entry qc_clustering \
   -profile apptainer \
   --input samplesheet.csv \
   --name my_study \
   --species human \
   --outdir results/qc_clustering
```

**Stage 2 — downstream analysis:**

```bash
nextflow run UKDRI/scdownstream -r dev_ukdri -entry downstream \
   -profile apptainer \
   --base_adata results/qc_clustering/my_study_finalized.h5ad \
   --name my_study \
   --species human \
   --outdir results/downstream
```

**Stage 3 — differential expression:**

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

> [!WARNING]
> Provide pipeline parameters on the command line or via `-params-file`. Custom config files passed
> with `-c` can supply any Nextflow configuration **except parameters**.

UK DRI users: see the
[UK DRI Informatics wiki](https://wiki.informatics.ukdri.ac.uk/en/Pipelines/nfcore_scdownstream) for
how to run the pipeline on the cluster.

## Changes and known limitations

This pipeline is **work in development**. The following are known and, for now, expected behaviours.
Several of them silently affect results, so please read before interpreting output.

1. **Doublets are flagged, not removed.** scrublet writes its scores and predictions into the object,
   but the doublet removal step is currently disabled. Filter on the doublet annotation yourself.
   `--doublet_detection_threshold` is not currently supported.
2. **Keep `scvi` in `--integration_methods`.** Everything after the merge is built on the scVI latent
   space; a selection that omits `scvi` produces no output and no error.
3. **Per-sample QC overrides in the samplesheet are ignored.** The per-sample `min_genes`,
   `max_mito_percentage`, `automatic_cell_filtering` and related columns are read into the sample
   metadata but the global `--min_*` / `--max_*` parameters are what actually reach the filtering
   step. Set filters globally.
4. **`--prep_cellxgene` is no longer supported.** Leave it at its default.
5. **`-profile test_offline` is no longer supported.** Use `-profile test` instead.
6. **Set `--species` explicitly.** It defaults to `human`, and mouse data analysed under the human
   default produces wrong enrichment and cell–cell communication results without any error.
7. **Five processes have no container registry fallback** and require locally built `.sif` images. We will publish these images to a container registry in the near future —
   see [Local container images](docs/usage.md#local-container-images).
8. **`--ortholog_hcop_directory` defaults to a UK DRI path** (`/nfsdata/genome/hcop/`). Off-site
   runs must override it.
9. **The legacy single-pass workflow is no longer supported** — always pass `-entry` (see the note
   above).
10. **`--unify_gene_symbols` is no longer supported.** HUGO-based gene symbol unification only
    applies to human data and is not reliable enough to recommend. Gene symbols are still harmonised
    across samples without it.
11. Several other inherited parameters are also not currently supported: `--skip_enrichment`,
    `--skip_liana`, `--skip_rankgenesgroups`, `--pseudobulk*`, `--cluster_per_label`,
    `--cluster_global`, and the `exclude_samples_col` / `exclude_samples_values` columns of the
    contrasts file.
12. MultiQC coverage is partial — the Quarto reports are the more complete view of a run.

## Documentation

- [Usage](docs/usage.md) — samplesheet format, all parameters, per-entry-point guidance
- [Output](docs/output.md) — the files each stage produces and how to read them
- [UK DRI Informatics wiki](https://wiki.informatics.ukdri.ac.uk/en/Pipelines/nfcore_scdownstream) —
  running the pipeline on the UK DRI cluster

Parameters are defined in [`nextflow_schema.json`](nextflow_schema.json); `nextflow run … --help`
lists them.

## Credits

This pipeline is derived from **nf-core/scdownstream**, originally written by
[Nico Trummer](https://github.com/nictru), and would not exist without the work of its authors and
contributors. Full attribution — original authors, upstream contributors, the nf-core framework and
the pipelines this work builds on — is in [`ACKNOWLEDGEMENTS.md`](ACKNOWLEDGEMENTS.md).

The fork is maintained by UK DRI Informatics.

## Contributions and support

Contributions are welcome via pull request against the `dev_ukdri` branch of
[UKDRI/scdownstream](https://github.com/UKDRI/scdownstream). For help, open an issue in this
repository or contact UK DRI Informatics.

## Citations

An extensive list of references for the tools used by the pipeline can be found in
[`CITATIONS.md`](CITATIONS.md).

This pipeline is no longer part of nf-core, but it was built on the nf-core framework and template,
which you may wish to cite:

> **The nf-core framework for community-curated bioinformatics pipelines.**
>
> Philip Ewels, Alexander Peltzer, Sven Fillinger, Harshil Patel, Johannes Alneberg, Andreas Wilm, Maxime Ulysse Garcia, Paolo Di Tommaso & Sven Nahnsen.
>
> _Nat Biotechnol._ 2020 Feb 13. doi: [10.1038/s41587-020-0439-x](https://dx.doi.org/10.1038/s41587-020-0439-x).
