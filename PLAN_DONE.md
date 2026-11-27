# UK DRI scdownstream fork: done log

This is the record of the work from [PLAN.md](PLAN.md) that is finished. The open work stays in
`PLAN.md`. Everything here is merged into `dev_ukdri` (up to PR #34, `bb5dd9f`) unless it says
otherwise.

> This document is now the **only** record of the container-URI history behind #4a (below), which
> the shipped docs deliberately no longer carry. (#4b, "the Dockerfiles are reconstructions", is
> resolved: the committed Dockerfiles are the recipes the published images were built from,
> pinned to the versions in the tested images.) Commit `0c5a914` removed both
> from `README.md` and `docs/usage.md` as development narrative. Do not prune them from here on
> the assumption they are written down elsewhere.

## Plan history

- **2026-09-01.** Partial in-place container remediation on the `readme` branch. The base image
  `gcfntnu/scanpy:1.11.4` was found to be public and enough for three processes, which removed
  two of the four hand-built images. Work on branch `containers` (`12102ea`, `0c5a914`).
- **2026-09-24.** Priority changed: a working, portable pipeline was needed before the new repo.
  Phase 0 became "make the fork usable by a general user in place".
- **2026-09-28.** All three stages (`qc_clustering` → `downstream` → `differential_genes`)
  complete on the UK DRI cluster with the public images, reports included.
- **2026-09-29.** Aim changed again: make a clean working version in this fork first, then port it
  into a new repo. `PLAN.md` was split into open work and this done log.
- **2026-10-21.** The cluster moved to Nextflow 26.04.4, whose strict parser rejected the fork
  (config errors, and `-entry` is no longer supported). Fixed in PR #19; stages are now selected
  with `--stage`. SCR_MM_1 and SCR_MM_2 then completed all three stages.

- **2026-11-06.** Custom model training merged: scANVI/scVI reference mapping in `qc_clustering`
  and `--stage train_models` (PR #30), and `--stage annotate` (PR #31).
- **2026-11-23.** PR #33 merged (`qc-max-thresholds-unset`, merge `8e32a86`): `--max_genes` /
  `--max_counts` unset by default (no upper limit) instead of `false`, schema `integer` (minimum 1),
  quality_control maps unset to `false` for SCANPY_FILTER; usage note "Upper thresholds". Fixes the
  `nf-core pipelines lint` abort on the mixed integer/boolean type. HPC run with defaults before
  the merge.
- **2026-11-23.** PR #34 merged (`hvg-seurat-v3`, merge `bb5dd9f`): Seurat v3 (vst) HVGs by
  default. `--hvg_flavor` / `--train_hvg_flavor` (`seurat_v3` | `seurat` | `pearson_residuals`) and
  `--hvg_span` / `--train_hvg_span` (0.3); PCA HVGs over all cells, integration and training HVGs
  per batch, on the raw counts (no rounding of ambient-corrected counts); clear errors for
  `n_hvgs 0` with seurat_v3, batches of < 2 cells and failed loess fits (fix: span 1); reports name
  and cite the method. HPC runs of qc_clustering and train_models passed; the trained model's
  balanced accuracy went from 0.988 (Pearson residuals) to 0.986 (seurat_v3).
- **2026-11-13.** PR #32 merged (`anndata-0.12-throughout`, merge `20bb8f0`; contains
  `marker-genes-presto-cosg`). All five stages passed on real data before the merge; the
  `seurat-ukdri-suite` modules (decontX, SoupX, scds, Seurat integration) were not run on real data.
  - **Marker genes:** presto statistics (Python re-implementation of `wilcoxauc`, equal to R presto
    1.1.0 to ~1e-14) and COSG; `--markers_method` (default `presto,cosg`), `--markers_rank_by
    auc|cosg`, `--markers_min_auc` (0.7), optional `--markers_min_cosg_score`;
    `--markers_thr_adj_pvalue` and `SCANPY_RANKGENESGROUPS` removed; marker JSON with `list_auc` /
    `list_cosg_score`.
  - **Why the containers changed:** COSG in the BioContainer (pandas 3.1 / anndata 0.13) wrote
    nullable string arrays that `SCANPY_ENRICH` (anndata 0.12.2) could not write back; then LIANA+
    (anndata 0.10) failed on the `null` elements that anndata >= 0.12 writes for `None`. Fix: one
    anndata generation (0.12.2) for all main h5ad steps: 23 modules on `gcfntnu/scanpy:1.11.4`, own
    images `nhecker/scanpy-ukdri-suite:1.11.4-1` (Python tools + Quarto, base packages as pip
    constraints) and `nhecker/seurat-ukdri-suite:5.5.1-1` (R 4.5.3, anndataR pinned 1.2.1 because
    1.2.2+ writes nullable strings; non-root). Dockerfiles in `containers/`.
  - **Removed:** `ADATA_TORDS` and the `.rds` outputs, SCimilarity with the unused pseudobulking path
    and the nf-core `untar` module. `ADATA_READRDS` adapted to anndata2ri 2.0.
    `utils_nfcore_scdownstream_pipeline` renamed to `utils`.
  - **Report:** marker column legend as a bullet list (AUC, COSG score, log2FC, pct, p, adj. p),
    marker/enrichment/LIANA+ tables sortable, numeric filter panel for the marker tables.

## Findings that are fixed

These are the original audit texts, with the fix below each.

### #4a. A Wave container URI was a placeholder that did not exist

*(verified against the registry, 2026-09-01)* Two were found:
`community.wave.seqera.io/library/pydeseq2_anndata:ed0f23c01ea418d9` and
`community.wave.seqera.io/library/decoupler_anndata:1f90c5f9ae5f4f8d`, both returning
`MANIFEST_UNKNOWN` (HTTP 404). All 46 other Wave references in the repo returned HTTP 200, as did
all five `docker.io/nicotru|saditya88` images and `quay.io/biocontainers/liana`, so this was not a
method artefact. Seqera Containers are documented as permanent (unlike Wave's default 1-week ECR
repo), so it was not expiry either. The commit that introduced them is titled *"Added pydeseq2
modules boiler plate"*. They were never built.

This was **worse than having no fallback**: the strings follow the correct Wave naming convention
and pass code review, and the modules ran fine on the cluster because
`params.singularity_cache_dir` shunted them to a local `.sif`. The failure appeared only off-site,
at pull time.

Fix:
- 2026-09-01: `DECOUPLER_PSEUDOBULK` and `FILTER_PSEUDOBULK` moved to `gcfntnu/scanpy:1.11.4`, a
  public Docker Hub image (miniconda3/bookworm, Python 3.11, scanpy 1.11.4) that already ships
  decoupler 2.1.1, anndata 0.12.2, pandas and pyyaml, the full dependency set of both templates.
  The pydeseq2 URI was replaced by the explicit marker `CONTAINER_REGISTRY_MISSING`.
- 2026-09-24: `DIFFERENTIAL_GENES_PER_CONTRAST` moved to `docker.io/nhecker/pydeseq2:0.1`, built
  from the committed Dockerfile (`5af9976`).
- The decoupler envs declared `decoupler>=1.6` while their templates call `dc.pp.pseudobulk` and
  `dc.pp.filter_samples`, which exist only in decoupler **2.x** (1.6 has `dc.get_pseudobulk` and
  no `dc.pp` namespace). Corrected to `>=2.0`. Treat it as evidence that the `>=` bounds were
  never resolved against the API the templates use (pinning them is still open, see `PLAN.md`).

### #4. Processes with no public container image

Originally four further processes had no public image: `scanpy/enrich`, `scanpy/report` ×2 and
`pydeseq2/report`. Three faults compounded:

- the fallback was an absolute site path (first `/data/nhecker/apptainer/images/…`, then
  `/nfsdata/apptainer/scanpy_1.11.4_coreinf_0.3.sif`);
- the ternary did not guard on `workflow.containerEngine`, so under `-profile docker` the
  container string resolved to a `.sif` path;
- several modules declared `conda "${moduleDir}/environment.yml"` with no such file, so
  `-profile conda` failed at parse.

Fix:
- 2026-09-01: `SCANPY_EXPORT_MARKERS` moved to `gcfntnu/scanpy:1.11.4`. `coreinf_0.1` and
  `coreinf_0.3` had no recoverable difference, so the remaining processes were collapsed onto
  `scanpy_1.11.4_coreinf_0.3.sif`.
- 2026-09-24: the report and enrichment processes moved to
  `docker.io/nhecker/scanpy-report:1.11.4-coreinf0.4`, built from the committed Dockerfile.
  Plain `container` directives, no site paths, no `params.singularity_cache_dir` (`5af9976`,
  `8d393ca`).
- `c482c2c`: Docker Hub images must name `docker.io/` explicitly, because `nextflow.config` sets
  the default registry to `quay.io`.
- `465e045`, `6cbf20f`: every module with a conda line now has an `environment.yml`, plus the
  `pydeseq2/report` conda line. (An earlier "no missing env files" check never matched anything,
  because of an unescaped `$` in the grep pattern. A fixed-string search over all 58 modules with a
  conda line finds no other gap.)
- `946024e`: Dockerfiles and conda envs pinned to the package versions in the published images.

### #6. `--ortholog_hcop_directory` shipped a site path as its schema default

It was `/nfsdata/genome/hcop/`. Fix (`31e3d50`): the default is `null`, and it is passed as `[]`
when unset. LIANA+ errors clearly only for non-human species. The docs point to
https://www.genenames.org/download/hcop/tsv/.

### #7. The default entry point was broken

`workflows/scdownstream.nf:153` read `params.clustering_resolution` (singular, never defined), so
a bare `nextflow run` passed `null` to Leiden. Fix (`3aa34e3`): the single-pass workflow was
removed; no `-entry` now warns and runs `qc_clustering`. (Its now-unused subworkflows and stale
pipeline tests are open, see `PLAN.md` #18.)

### #13 (part). `manifest.defaultBranch = 'master'`

No `master` branch exists, so `nextflow run UKDRI/scdownstream` without `-r` failed. Fix
(`bdad8e6`): `defaultBranch = 'dev_ukdri'`, which is also GitHub's default branch. The rest of #13
(pipeline identity) is new-repo work.

### #21. `downstream` and `differential_genes` required `--input`

They never read the samplesheet, so users passed a dummy file. Now not required; confirmed in
the user's runs.

### #1. Mitochondrial filtering ran on a metric from another module (closed 2026-10-22)

*Original:* `--mito_genes` was accepted and never used, and `SCANPY_FILTER` filtered on
`obs["pct_counts_mt"]` (and the other QC columns) computed by QC_RAW, before ambient correction, so
skipping or reordering QC_RAW would give a `KeyError`.

*Closed:* `--mito_genes` was removed and moved to the backlog (PR #25). The
`KeyError` case can no longer happen: QC_RAW runs unconditionally in `QUALITY_CONTROL`, and there is
no option to skip QC (`--qc_only` was removed in PR #22). Computing the QC metrics on the counts
before ambient correction is kept as a design choice: with the default `--ambient_correction none`
they are the same counts, and with correction on, filtering on raw metrics is a common choice.
Remaining edge: the integration-time `SCANPY_FILTER` also reads those obs columns, so a
`--base_adata` not produced by this pipeline must carry `total_counts`, `n_genes_by_counts` and
`pct_counts_mt`. Mito QC without gene symbols is tracked as #48.

### #2. Per-sample QC overrides were computed and then discarded

*Original:* `quality_control/main.nf` built a `multiMap` of per-sample thresholds from `meta`, then
called `SCANPY_FILTER` with `params.*`; the samplesheet QC columns were dead. Decision D2.

*Fix (PR #26, `2ab7aab`):* each threshold is the sample's samplesheet value if set, otherwise the
global parameter (`sampleOrGlobal()` in `utils`), passed to
`SCANPY_FILTER` as per-sample channels from the same `multiMap`. Automatic filtering (per sample or
global) still replaces the fixed thresholds, with `min_genes` / `min_counts_cell` as floors. The
schema defaults of the QC columns were removed (they hid "unset"; `automatic_cell_filtering`
defaulted to `true`), `max_genes` / `max_counts` columns added (minimum 1), and the old `?:` logic,
which dropped explicit `0` / `false` and read `min_counts` instead of `min_counts_cell`, replaced.
Checked with nf-schema 2.7.3: an empty samplesheet cell arrives as `[]` and falls back to the global.

### #5. `-profile test_offline` included a file that did not exist

Removed in PR #19 (`conf/test_offline.config` was never committed). The strict parser checks
include paths even for unselected profiles, so it broke every run under Nextflow 26.04.

### #8. `--prep_cellxgene true` crashed `downstream` and `qc_clustering`

Option and `ADATA_PREPCELLXGENE` removed in PR #22. Besides the wiring bug, the template would
have log1p-transformed the already log-normalised `X` and stripped the `X_` prefix cellxgene
needs. The published h5ad opens in cellxgene directly.

### #11. `conf/test_full.config` set the undeclared `celldex_reference_label`

Removed with SingleR in PR #20.

### #14 (part). Duplicate parameters

`celltypist_model` declared once (`null`) and `celldex_reference` removed, PR #20. `outdir` and
`email` are still duplicated (open in `PLAN.md`).

### #17. `exclude_samples_col` / `exclude_samples_values` were a fully plumbed no-op

*Original:* validated by `schema_contrasts.json`, parsed in `main.nf`, carried through
`workflows/differential_genes.nf`, never read by `DIFFERENTIAL_GENES_PER_CONTRAST`. Decision D3.

*Fix (PR #24, `cc9981c`):* implemented with nf-core/differentialabundance semantics. For that
contrast only, `deseq2_diffgenes_per_contrast.py` drops every pseudobulk sample whose value in
`exclude_samples_col` is one of `exclude_samples_values`, before the counts, the design and the
`--diffgenes_min_samples` check. No workflow change was needed: decoupler's pseudobulk (2.1.1,
checked) keeps every obs column that is constant within a profile, so sample-level columns reach
the DE step; adding the column to the pseudobulk grouping columns was rejected because a
non-sample-level column would split the profiles of every contrast. A missing column stops with an
error explaining it must be sample-level; unmatched values and half-set pairs are warnings; a skip
caused by the exclusion says so. Tested by rendering the template with pydeseq2 0.5.4 / anndata
0.12.2: excluding two donors gives the same results as a pre-subset input.

### #15. `skip_enrichment` was in `nextflow.config` but not the schema

Removed in PR #25 (nothing read it); an enrichment skip option is in the `PLAN.md` backlog.

### #18 (part). Dead code

`subworkflows/local/{cluster,per_group}` deleted in PR #19, the unused `FINALIZE` subworkflow in
PR #22. The rest of #18 is open.

### #24. Ambient RNA correction ran by default on already-corrected input

Default `--ambient_correction none` (D4), PR #17 (`eb8c069`, docs `eda74f3`).

### #26 (part). Stubs with an empty `versions.yml`

`ADATA_PUBLISH`'s stub writes a real `versions.yml` (PR #21). 42 stubs remain (open).

### #27. Nextflow 26.04's strict parser rejected the config

PR #19: `test_offline` removed, the `validation.help` / `summary` text that read `manifest.*` and
`validation.help.*` dropped, `ANNDATA_GETSIZE` no longer declares `prefix` twice, and nf-schema
bumped 2.4.2 → 2.7.3 (its `validation` options are unknown to the strict parser in 2.4.2), which
raises the minimum Nextflow to 25.10.0. CI tests 26.04.4. `nextflow lint` reports 0 errors.

## Feature: clearer differential_genes report (done, PR #27)

The planned feature as recorded in `PLAN.md`, kept for reference:


*Recorded 2026-10-22; implemented 2026-10-22 on branch `de-report-improvements` (all four changes plus the overview table), rendered and checked on a synthetic data set.* Report: `modules/local/pydeseq2/report/templates/scdownstream_differential_genes_report.qmd`,
called by `PYDESEQ2_GENERATE_REPORT` (`workflows/differential_genes.nf:118`).

#### Findings (open)

**44. The contrasts aren't stated.** The report only sees contrast ids (third column of the DE
manifest, `MANIFEST`). The contrasts file (`--diffgenes_contrasts`: id, variable, reference, target,
blocking, exclude_samples_col/values) is not an input of `PYDESEQ2_GENERATE_REPORT`, so a reader
can't tell what `AD_vs_ctrl` compares or in which direction the fold change goes.

**45. "Group" doesn't say which column.** Sections are headed `Group: <value>` (`.qmd:465`). The
column name is available (first entry of `GROUP_COLUMNS`, i.e. `--diffgenes_group_col`) but isn't
shown, so "Group: Microglia" doesn't tell the reader it is e.g. `celltypist:majority_voting`.

**46. Volcano plots with no significant genes take the same space** as informative ones
(`.qmd:489-501`). They are drawn and shown in full whenever results exist.

**47. Up/down counts are hard to find.** The number of significant genes is only in the volcano
figcaption ("n of m tested gene(s) pass…", `.qmd:501`), and is not split into up and down.

#### Proposed changes

- **Contrasts table at the top** (before the pseudobulk summary). Columns: contrast id, variable,
  target vs reference, with the direction of log2FC spelled out ("positive = higher in target"),
  blocking variables and sample exclusions. This needs the contrasts file as a new report input.
  - **Interface change** (decide first, per CLAUDE.md): add `path(contrasts)` to
    `PYDESEQ2_GENERATE_REPORT`, wired from `params.diffgenes_contrasts` in
    `workflows/differential_genes.nf`, plus a `CONTRASTS_FILE` qmd parameter. nextflow-dsl2
    handles the wiring; pipeline-developer handles the module and qmd.
- **Name the group column everywhere:** state it in the intro ("Groups are the values of
  `<column>`"), and head the sections `<column> = <value>`, not `Group: <value>`.
- **Overview table at the start of each group section:** one row per contrast with significant
  **up**, significant **down** and tested genes, using the same `padj`/`|log2FC|` thresholds as the
  volcano plots. Optionally also one overall group × contrast table near the top of the report.
- **Collapse empty volcano plots:** when a contrast has no significant genes, wrap the plot in
  `<details><summary>Volcano plot (no significant genes) — click to show</summary>…</details>`. Plots
  with hits stay expanded.
- Update `docs/output.md` (report description) and the CHANGELOG.

#### Verification

Render the qmd with a small pseudobulk object, manifest and contrasts file (two groups, two
contrasts, one of them without hits). Check the contrasts table, the column name in the headings,
the up/down table and that only the empty volcano is collapsed.


## Feature: custom model training and cell type annotation (done, PRs #30, #31)

*Planned 2026-09-29; Phase 1 implemented 2026-10-28 on `scanvi-query`, Phase 2 on `train-models`
(stacked); `--stage annotate` added 2026-11-05 on `annotate-stage`. All merged 2026-11-06. Open
findings and later items stay in [PLAN.md](PLAN.md).*

### Decisions

- **A new stage in this repo, `--stage train_models`,** not a separate pipeline. A trained
  model must be used with the same tool versions (scvi-tools 1.3.3, celltypist 1.6.3), containers
  and gene preprocessing it was trained with. Loading, `UNIFY`, HVGs and the scvi-tools modules
  are shared with `qc_clustering`. A separate repo only pays off if training becomes its own
  product (large atlases on GPU, its own release cycle, used by other pipelines).
- **Training input:** an external annotated reference h5ad, or this pipeline's own stage output
  after the user has curated labels (counts in `layers["counts"]`).
- **Training gets its own preparation step, and keeps labels exactly as given** (decided
  2026-10-06), as part of `ADATA_TRAINSPLIT` or a small `ADATA_TRAINPREP` before it: raw counts
  from `--train_counts_layer` into `X`, `--train_label_col` → `label` with only
  `--train_unknown_label` and NaN mapped to `unknown`, `--train_batch_col` → `batch`, gene symbols
  from `--train_symbol_col`. `ADATA_UNIFY` (used only by `qc_clustering`,
  `qc_clustering.nf:133`) is **not** reused: its label normalisation
  makes distinct labels collide (#37), and its gene/batch harmonisation is meant for merging
  samplesheet samples.
- **All new parameters are prefixed `train_*`.** The per-stage schemas (A3) come later.
- **Two phases, each its own branch and PR:** `scanvi-query` (Phase 1), then `train-models`
  (Phase 2), both off `dev_ukdri`.
- **`--base_adata` together with scANVI** is blocked with a clear "not supported yet" error for
  now (see #31).
- **Prediction columns:** `scanvi:<model>` and `scanvi:<model>:conf`, the same convention as
  `celltypist:<model>` / `celltypist:<model>:conf`.
- Implemented through the pipeline-developer and nextflow-dsl2 agents (`CLAUDE.md`).

### Phase 1: scANVI models as input in `qc_clustering` (branch `scanvi-query`)

Target: `--integration_methods scvi,scanvi --scanvi_model <model>` maps new cells onto the
reference (no `--base_adata`) and adds `scanvi:<model>` / `scanvi:<model>:conf`, used by the
reports exactly like CellTypist columns.

- **pipeline-developer:**
  - `scvitools/scanvi` and `scvi`:
    - use an explicit mode input (`scratch` / `from_scvi` / `query`) instead of guessing from
      `meta2.id` (#28);
    - use `layer="counts"` everywhere, creating it from `X` when missing (#32);
    - for queries, set labels to the model's unlabelled category (read from the model registry),
      keep the input labels as `label_input`, and use `weight_decay=0.0` with
      `task.ext.query_max_epochs` (#33, #34);
    - `predict(soft=True)`: argmax gives `scanvi:<model>`, max probability gives
      `scanvi:<model>:conf`; write `X_scanvi` (#35, #38);
    - move seed, precision and early stopping into named variables (repo rule).
  - Reports: select label columns by a list of prefixes (`celltypist:`, `scanvi:`), and handle
    non-integer groups (#39).
- **nextflow-dsl2:**
  - `integrate.nf` / `combine.nf` / `qc_clustering.nf`:
    - pass the mode;
    - give scANVI a query input with all genes and counts whenever a scANVI model is given
      (#36);
    - stage models with `file(..., checkIfExists: true)`;
    - route scANVI obs/obsm into the object that reaches `FINALIZE_H5AD_SIMPLE` (#29).
  - New param `--scvi_query_max_epochs` (default 100).
  - Checks:
    - `--scanvi_model` requires `scanvi` in the methods;
    - `scanvi` requires `scvi` (#30);
    - `--base_adata` with `scanvi` stops with "not supported yet" (#31).
- **Docs:** rewrite "Reference mapping" to match the new behaviour.

### Before Phase 2

Checked 2026-10-06 against `dev_ukdri` `ca09553`:
- **Phase 1 complete.** Its blockers still hold (#28, #29, #32, #33, #35, #39), and Phase 2
  evaluates scANVI by mapping the test cells as a query.
- **Needed in Phase 2:** `SCANPY_HVGS` gets a selectable flavour, so the reference model can use
  `seurat_v3` on counts with a batch key (`batch_key` is already supported).
- **Useful:** #26, so that stubs write a valid `versions.yml` and `ADATA_EXTEND` gets a stub; then
  `--stage train_models -stub` can test the wiring locally.
- **Not needed first:** #37 (training does not use `UNIFY`), the doublet clean-up, the per-stage
  schemas (`train_*` prefixes instead), and `--base_adata` mapping (blocked).

### Phase 2: `--stage train_models` with evaluation and report (branch `train-models`)

```
--train_adata ─► ADATA_TRAINSPLIT (prepare + split) ─► all / train / test (label hidden) + label_counts.tsv
 evaluate:  CellTypist: CELLTYPIST_TRAIN(train) ─► CELLTYPES_CELLTYPIST(test, .pkl) ─┐
            scANVI: HVGS(seurat_v3, batch) ─► SCVI(train) ─► SCANVI(from_scvi, train) │
                    ─► SCANVI(query: test) ──────────────────────────────────────────┤
                                                           MODEL_EVALUATE per method ◄┘
 final:     CELLTYPIST_TRAIN(all); HVGS ─► SCVI(all) ─► SCANVI(all)
 report:    SCANPY_GENERATE_REPORT_TRAINING
```

**Evaluate, then refit on all cells.** Performance is estimated on a held-out test set, and the
final model is trained on all cells with exactly the same settings. The test metrics then describe
the training *procedure*; they are slightly pessimistic for the final model, which sees more data.
- **Small labels:** labels with fewer than `--train_min_cells_per_label` cells are dropped before
  both the split and the final fit, and listed in the report. Otherwise the final model would know
  classes the evaluation never tested.
- **Stratified split:** `--train_test_fraction` of each label goes to the test set, with a fixed
  seed held in a named variable.
- **Cell-level splits are optimistic for new donors,** because test cells share samples with the
  training cells. `--train_test_split_by <obs column>` (e.g. donor) holds out whole groups instead,
  via `StratifiedGroupKFold`, and warns for labels that end up with no test cells.
- **scANVI is evaluated the way it is used:** the reference is trained on the training split, and
  the test cells are mapped as a query with their labels hidden (Phase 1 code).
- **Unlabelled (`unknown`) cells** are never in the test set. CellTypist excludes them; scANVI
  keeps them as unlabelled cells for semi-supervised training.
- **Compute cost:** every model is trained twice, which matters for scANVI on large atlases, so
  there are switches to skip the refit or the evaluation.

| Parameter | Default | Meaning |
|---|---|---|
| `--train_adata` | — (required) | Annotated reference h5ad |
| `--train_label_col` | — (required) | `obs` column with the cell type labels |
| `--train_unknown_label` | `unknown` | Label value treated as unlabelled |
| `--train_batch_col` | `batch` | Batch column; missing → one batch |
| `--train_counts_layer` | `X` | Where the raw counts are; `counts` for this pipeline's outputs |
| `--train_symbol_col` | `index` | `var` column with gene symbols |
| `--train_methods` | `celltypist,scanvi` | Models to train |
| `--train_name` | `--name` or `reference` | Model file names |
| `--train_min_cells_per_label` | `50` | Labels with fewer cells are dropped everywhere |
| `--train_test_fraction` | `0.25` | Share of each label held out for testing |
| `--train_test_split_by` | none | `obs` column to hold out whole groups instead of cells |
| `--train_evaluate` | `true` | Run the split and the evaluation |
| `--train_final_on_all` | `true` | Refit the final model on all cells; `false` keeps the training-split model |
| `--train_hvgs` | `5000` | HVGs for the scVI/scANVI reference model |
| `--train_celltypist_feature_selection` | `true` | CellTypist two-pass feature selection |
| `--train_celltypist_use_sgd` | `false` | SGD / mini-batch training for large references |

scVI/scANVI architecture and epochs reuse `--scvi_*`, because a query inherits the architecture
from the saved model.

**New modules (pipeline-developer):**
- `ADATA_TRAINSPLIT`: drops small labels and makes the split. It emits `all.h5ad`, `train.h5ad`,
  `test.h5ad` (true label moved to `label_true`, `label` set to `unknown`) and `label_counts.tsv`
  (label, n_total, n_train, n_test, dropped).
- `CELLTYPIST_TRAIN`: normalises with `normalize_total(1e4)` + `log1p` (matching
  `celltypist.py`) and calls `celltypist.train` (arguments verified against 1.6.3).
- `MODEL_EVALUATE`: overall accuracy, balanced accuracy and macro F1; precision, recall, F1 and
  support per label; confusion matrix. Written as TSVs.
- `CELLTYPES_CELLTYPIST` accepts local `.pkl` files as well as model names.
- `SCANPY_GENERATE_REPORT_TRAINING` + `scdownstream_training_report.qmd`:
  - label counts (train, test, dropped);
  - per method: metrics, a per-label table and a confusion matrix normalised per true label;
  - a note on what the split means.

**Wiring (nextflow-dsl2):**
- `workflows/train_models.nf`;
- the `train_models` entry in `main.nf` (no `--input`, completion summary, error checks);
- the `train_*` params in a "Model training options" schema group;
- `--celltypist_model` accepting `.pkl` paths mixed with names, split in `celltype_assignment`.

**Outputs:** under `${outdir}/models/`:
- `celltypist/<name>.pkl`;
- `scanvi/<name>/model.pt` and the scVI model;
- `<name>_reference.h5ad`, with the scANVI latent and labels;
- `evaluation/<method>/*.tsv`;
- the training report.

**Using the models in `qc_clustering`:** `--celltypist_model my_ref.pkl` (mixed with built-in
names as needed), and `--scanvi_model models/scanvi/<name>/model.pt` with
`--integration_methods scvi,scanvi`.

### `--stage annotate` (PR #31)

- Adds CellTypist / scANVI predictions from pre-trained models to a processed h5ad
  (`--base_adata`) without re-running QC, integration and clustering; models via
  `--celltypist_model` / `--scanvi_model`; `--annotate_batch_col`, `--annotate_symbol_col`, and a
  cell subset (`--annotate_subset_col` / `--annotate_subset_values`, output = subset only); small
  HTML report; predictions table as TSV (no `.pkl` hand-offs).
- Follow-up `eeeb6b6` after a crash on an input without a `counts` layer: raw counts from
  `layers["counts"]`, else from `X`; no log-normalisation or integer check (raw counts are a
  documented requirement; ambient correction gives non-integer counts); batches from
  `--annotate_batch_col`, else one batch; only this run's prediction columns are dropped from the
  working copy and copied back; CellTypist takes the raw counts from `layers["counts"]` itself.

### Verification

- **Phase 1, locally** (scratch venv with scvi-tools 1.3.3 on CPU):
  - train a tiny scVI → scANVI reference, save it, and run the rendered `scanvi.py` in query mode
    on a second object with unseen and `unknown` labels;
  - check that the model is loaded rather than retrained, the labels are hidden, and the
    `scanvi:*` columns and `X_scanvi` are written;
  - `-preview` of `qc_clustering` with and without `--scanvi_model`, plus the new checks;
  - run the report cells on an object with both `celltypist:` and `scanvi:` columns.
- **Phase 2, locally:**
  - `ADATA_TRAINSPLIT` on a synthetic object with a small label, `unknown` cells and a donor
    column; label names stay unchanged ("CD4+ T cells" and "CD4- T cells" remain two classes);
  - CellTypist train, then predict with the `.pkl` path;
  - `MODEL_EVALUATE` metrics checked against sklearn by hand;
  - the report cells;
  - `-preview` of every switch combination.
- **Both phases, on the cluster:** train on a stage-1 test output (`--train_counts_layer counts`,
  a CellTypist column as stand-in labels, a lower `--train_min_cells_per_label` for the tiny test
  data), check the training report, then use both models in `qc_clustering`.

## Other fixes found while testing

Testing ran the test profile end to end on the UK DRI cluster (Apptainer), from 2026-09-24:

```bash
nextflow run https://github.com/UKDRI/scdownstream -r dev_ukdri -latest \
  -entry qc_clustering -profile test,apptainer --outdir results_qc -resume
```

| Step | Failure | Fix | Where |
|---|---|---|---|
| `SCANPY_SCRUBLET` | h5ad write `TypeError` on `predicted_doublet`: a cell scrublet can't score gets `NaN`, and the write came before the `bool` cast (UK DRI edit in `9e7be03`) | Cast before writing | PR #7 (`0156d41`) |
| `CELLTYPES_CELLTYPIST` | Empty (0, 0) AnnData: the default `min_genes 200` / `min_cells 5` removed every cell of the shallow nf-core test sample SRR28679758 (median ~27 genes per cell) | `min_genes 10`, `min_cells 3` in `conf/test.config` only | PR #8 (`62f9a07`) |
| Container pull | Docker Hub images pulled as `quay.io/...`, because `nextflow.config` sets the default registry to `quay.io` | Prefix all Docker Hub images with `docker.io/` | PR #9 (`c482c2c`) |
| `SCANPY_HVGS` | Missing `*_hvg.pkl`: not written when the object has ≤ `n_hvgs` genes | Always write it; mark every gene highly variable in that case | PR #9 (`f8a3135`) |
| QC/clustering report | `KeyError: 'hvg'` in `sc.pl.highly_variable_genes`: the fallback above skipped the HVG statistics | Always run `sc.pp.highly_variable_genes` | PR #9 (`79fd8de`) |
| Analysis report (LIANA+ dotplot) | `['9'] not found in liana_res['target']`: LIANA+ drops clusters with < 5 cells, and the dotplot requested every cluster | Only plot clusters present in `liana_res`; note the skipped ones | PR #11 (`1184fab`) |

**Result (2026-09-28):** `qc_clustering` (test profile), then `downstream` from its output, then
`differential_genes` (with the user's own contrasts) all complete on the UK DRI cluster with
Apptainer and the public images, reports included.

## Branches and PRs

**`containers`** (PR #6 `7d8952e`, PR #9 `2e0d586`, branch deleted): the container work above,
`12102ea`, `0c5a914`, `5af9976`, `8d393ca`, `c482c2c`, `f8a3135`, `79fd8de`.

**`general-user-fixes`** (off `2e0d586`, 2026-09-25; PR #11 `df981b0`, branch deleted). Reviewed
by the nextflow-dsl2 and pipeline-developer agents, no blocking issues.
- `31e3d50`: optional HCOP directory (#6).
- `bdad8e6`: `defaultBranch = 'dev_ukdri'` (#13).
- `465e045`: the missing conda `environment.yml` files, plus the `pydeseq2/report` conda line (#4).
- `3aa34e3`: broken default workflow removed; no `-entry` runs `qc_clustering` (#7).
- `d28a17e`: stubs and versions for `scanpy/report`, `adata/add_obs_obsm` and
  `liana/rankaggregate` (this restored the report modules' `versions.yml` output, part of #20).
- `946024e`: Dockerfiles and conda envs pinned to the published image versions. No rebuild needed.
- `59ea8db`, `e8a7d3e`: the README says the images are pulled from Docker Hub, and has a
  Prerequisites section (Nextflow, Apptainer, x86_64 Linux, internet; HCOP for non-human species).
- `1184fab`: LIANA+ dotplot fix.

**`default-scvi-umap`** (off `df981b0`, 2026-09-28; PR #12 `681d7ec`, branch deleted).
- `8f02632`: the reports plot the scVI UMAP (`X_umap_scvi`) by default, with the PCA UMAP as the
  secondary view, since clustering uses the scVI graph. Not changed, because they draw no plots in
  the current stages: the Leiden module plot (`plot_umap=false`) and `ENTROPY`
  (`ext.plot_basis = 'umap'`, never called).
- `9422f49`: `--umap_for_plots` (default `X_umap_scvi`), passed to both reports with `-P`. Missing
  keys fall back to `X_umap_scvi` → `X_umap_pca` → `X_umap` with a note in the report. The PCA UMAP
  is stored as **`X_umap_pca`**. `X_umap` holds the same embedding, kept as a copy (via
  `ext.alias_key` on the qc_clustering `SCANPY_UMAP` call) only because functions such as
  `sc.pl.umap` require it. `docs/output.md` lists the obsm keys: `X_scvi`/`X_emb`, `X_umap_scvi`,
  `X_pca`, `X_umap_pca`, `X_umap`.

**`output-naming-gzip`** (2026-09-28; PR #13 `67bca22`, branch deleted).
1. `6bea79e`: stage outputs are named after `--name` and the stage (fallback `scdownstream`):
   `<name>_qc_clustering.h5ad`/`.rds`, `<name>_downstream.h5ad`/`.rds`,
   `<name>_downstream_markers.json.gz`, and the reports `<name>_qc_clustering_report.html` /
   `<name>_downstream_report.html`. Before, stage 1 was `integrated_scvi_finalized.*`, because
   `integrate.nf:91` resets `meta.id`, so `--name` never reached it. The prefixes use
   `params.name` directly in `conf/modules.config`. Underscores rather than dots, because
   Nextflow's `simpleName` cuts at the first dot.
2. `bf4e7f5`: `SCANPY_UMAP` had no publish rule and was **always** published to
   `${outdir}/scanpy`, with the scVI and PCA calls (same `meta.id` `integrated_scvi`) overwriting
   each other. `SCANPY_NEIGHBORS` clashed the same way with `--save_intermediates`. Fixed with
   per-call prefixes (`_umap_SCVI`/`_umap_PCA`, `_neighbors_SCVI`/`_neighbors_PCA`) and a
   `SCANPY_UMAP` publish block (`${outdir}/scanpy/${meta.id}/`, only with `--save_intermediates`).
3. `2d841a1`: `ADATA_PUBLISH` rewrites the final h5ad with `write_h5ad(compression="gzip")` in
   `docker.io/gcfntnu/scanpy:1.11.4`, instead of `mv`. All five h5ad-writing images (gcfntnu,
   anndata_pyyaml, pyyaml_scanpy, scanpy-report, pydeseq2) have h5py 3.13.0 / HDF5 1.14.3 with the
   deflate filter. On nf-core SRR28679756: **5.2 MB → 1.8 MB (~65 % smaller)**, identical content,
   backed reads work, 0.1 s write. Only the published file is compressed. It guards against
   writing through the input symlink.
4. `6cbf20f`: the last missing conda env (`scanpy/export_markers`).

**`add_metadata`** (2026-09-29; PR #14 `ff50880`, branch deleted).
- `ADATA_ADDSAMPLE` sets `obs["sample"]` right after loading; a differing existing column is
  kept as `sample_original`. `ADATA_ADDMETADATA` adds the columns of a per-sample TSV
  (`--metadata`, `--metadata_sample_col`) before QC. It stops on a missing sample, duplicate ids,
  duplicate or empty header names, a column that already exists in `obs`, and the reserved names
  `batch`, `label` and `sample_original`. Bool columns are stored as string categories.
- `ADATA_MERGE` always keeps `sample` and the metadata columns (`force_obs_cols`), also with
  `--base_adata`.
- Fixed the `ADATA_UNIFY` sample check: `Series.equals()` against a string was always False.
- `--umap_color_by` / `--umap_color_by_embeddings`: both reports get a "Requested UMAP plots"
  section, one block per UMAP (PCA = before, scVI = after integration), one panel per column.

**`pca-n-pcs`** (2026-09-29; PR #15 `d8267dd`, branch deleted).
- `--pca_n_comps` (PCs computed, default 50) and `--neighbors_n_pcs` (PCs used for the PCA graph
  and the PCA UMAP, default all), passed via `task.ext` to the `QC_CLUSTER` calls only. `pca.py`
  caps `n_comps` at `min(n_obs, genes used for PCA) - 1` (the highly variable genes when set),
  as scanpy did itself when `n_comps` was not passed; startup check `neighbors_n_pcs <=
  pca_n_comps`. Elbow plot in the QC/clustering report.
- `CLAUDE.md` and the `.claude/agents/` definitions committed; new pipeline-developer rule
  against hard-coded parameters in function calls.

**`keep-outliers`** (2026-10-01; PR #16 `ca09553`, branch deleted; tested on GSE307000).
- `--filtering_keep_outliers`: cells failing the QC thresholds are kept and marked in a bool
  `obs["outlier"]` column. Only the `QUALITY_CONTROL` call of `SCANPY_FILTER` gets it
  (`withName: '.*:QUALITY_CONTROL:SCANPY_FILTER'`); the `INTEGRATE` call still drops zero-count
  cells before scVI. `filter.py` builds one mask with the old criteria (remove mode unchanged),
  and the early filtering pass in fixed-threshold mode is gone, so the histograms show all cells.
  The QC/clustering report gets a "QC outliers" part (UMAP, per-sample table, QC-metric violins).
- `--doublet_removal` (default off): removes cells called doublets by at least
  `--doublet_detection_threshold` methods right after detection, via `DOUBLET_REMOVAL`, plus a
  "dedoubleted" row in the cell-count table. `DOUBLET_REMOVAL` aligns predictions to the cells,
  counts missing/NA predictions as not doublets, skips the upset plot when no intersection reaches
  its minimum size (upsetplot crashed on an empty plot), names methods by file, and has a pinned
  `environment.yml` (today's unpinned resolve broke upsetplot 0.9.0). Prefix `_dedoubleted`.
- Resolved decision D5: scrublet's unscorable cells (`NaN`) count as not doublets.
- Fixed: `--automatic_cell_filtering` removed every cell as soon as one cell had zero counts
  (`pct_counts_mt` NaN made the thresholds NaN); thresholds now ignore NaN. Also the lower-bound
  clamp typo in `get_thresholds()`.
- Fixed: `ADATA_MERGE` turned bool `obs` columns into "True"/"False" strings; they stay bool when
  every object has them. `DOUBLET_DETECTION` reported only the last method's versions.
- `outlier` added to the reserved `--metadata` column names.

**#4b resolved (2026-10-06).** The Dockerfiles in `modules/local/scanpy/report/` and
`modules/local/pydeseq2/differential_genes/` were thought to be reconstructions. The user
confirmed they match the recipes used to build `docker.io/nhecker/scanpy-report:1.11.4-coreinf0.4`
and `docker.io/nhecker/pydeseq2:0.1`, apart from the version pins added in `946024e` (read from
the published images; the original recipes installed mostly unpinned versions, and their
unquoted `pydeseq2>=0.4` was a shell redirection). The images are tested, so the planned
image-parity check against the old cluster `.sif` files is dropped.

**`ambient-default-none`** (2026-10-07; PR #17 `e101d8c`, branch deleted). #24 / D4.

**`diffgenes-markers`** (2026-10-20; PR #18 `5a3bd54`, branch deleted). A marker-gene JSON for
`differential_genes`, `<name>_differential_genes.json.gz`, with set ids built with a configurable
`--markers_sep`; the default `--` followed on `dev_ukdri` (`6ae8cd9`).

**`fix/nf26-strict-parser`** (2026-10-21; PR #19 `9652d88`, branch deleted). #5, #27, `--stage`
replacing `-entry`, a nested-closure bug in `SINGLER` that only showed at runtime under the strict
parser, unused `CLUSTER` / `PER_GROUP` subworkflows deleted.

**`remove-singler`** (2026-10-21; PR #20 `325fa4e`, branch deleted). SingleR / celldex annotation
removed (no run used it; annotation is CellTypist, scANVI planned). #11, #14 (part).

**`publish-sparse-gzip`** (2026-10-22; PR #21 `c85fd44`, branch deleted). `ADATA_PUBLISH` converts
`X`, layers and `raw.X` to CSR before its gzip write, so `<name>_qc_clustering.h5ad` and
`<name>_downstream.h5ad` are always sparse and gzip-compressed. Stub fixed (#26 part).

**`fix-finalize-qconly-unify`** (2026-10-22; PR #22 `d02adec`, branch deleted). `--prep_cellxgene`
(#8) and `--qc_only` removed (the latter published every per-sample object under one name), the
unused `FINALIZE` subworkflow deleted, and `ADATA_UNIFY`'s duplicate-gene aggregation vectorised:
`--duplicate_var_resolution max` crashed on CSR input.

**`fix-empty-cells-filtering`** (2026-10-22; PR #23 `bc4a95f`, branch deleted). GSM5223974: the N-MAD lower bounds
clamped to 0 and a 1-UMI cell became empty after gene filtering. `--min_genes` / `--min_counts`
are now floors in automatic mode, and cells with no counts after gene filtering are always
removed (also with `--filtering_keep_outliers`).

**`contrast-exclude-samples`** (2026-10-22; PR #24 `9c68f80`, branch deleted). #17 / D3.

**`remove-mito-genes`** (2026-10-22; PR #25 `688c686`, branch deleted). `--mito_genes` and
`--skip_enrichment` removed (both unused; backlog), `SCANPY_FILTER` lost its `mito_genes` input, its
module tests now match the signature. #15, and #1 closed.

**`per-sample-qc`** (2026-10-22; PR #26 `af45a50`, branch deleted). #2 / D2, plus `--symbol_col`
removed (unused; gene symbols are the per-sample samplesheet `symbol_col`) and a "Gene symbols"
section in `docs/usage.md` documenting `symbol_col` and the limits of `symbol_col = none`.

**`de-report-improvements`** (2026-10-22; PR #27 `4c21ae8`, branch deleted). Findings #44-#47:
the differential_genes report gets a Contrasts table (direction of log2FC in words, blocking,
exclusions; new `PYDESEQ2_GENERATE_REPORT` input written by the workflow from the validated
contrasts), a Significant genes overview, `<group column> = <value>` headings, a per-group up / down
/ tested table before the plots, collapsed volcano plots without hits, and one significance rule
throughout. Rendered with Quarto 1.10 on a synthetic data set (real PyDESeq2 results); the render
caught an `inline_table_html` bug (pandas' indentation became an indented code block in pandoc),
fixed before merge.

**`fix-liana-dotplot-targets`** (2026-10-23; PR #28 `06dcab5`, branch deleted). SCR_HS_1's
downstream report crashed with "['40'] not found in liana_res['target']": liana 1.6.1's
`li.pl.dotplot` filters by source first and then requires every requested target among those rows,
but the report passed one global target list. Targets are now taken per source cluster, and a
failing dotplot becomes a note (the cluster's table is still shown). Tested by running the report's
LIANA cells with liana 1.6.1 on a synthetic `liana_res`.

**`report-top-celltype-labels`** (2026-10-23; PR #29 `a619bfa`, branch deleted).
`--report_max_celltype_labels` (default 50): the "Predicted labels" UMAPs and count tables of the
QC/clustering and downstream reports show the top N CellTypist labels by number of cells, the rest
grouped as grey `other`; models with ≤ N labels unchanged; the per-cluster label summary
(top 5 labels ≥ 5 % per cluster) stays unlimited. Tested on synthetic CellTypist output with
scanpy 1.11.4.

**`scanvi-query`** (`e23a0ef`) and **`train-models`** (`8e95f4d`, `1f7b89a`, `5e19995`; PR #30
`ed67252`, branches deleted 2026-11-06): scANVI/scVI reference mapping in `qc_clustering`,
`--stage train_models`, the params-file advice, and the QC report fix for ≤ 3 clusterings
(`5e19995`, cherry-picked from **`fix-qc-report-cluster-grid`** `797f3c4`, whose own PR was not
needed; branch deleted).

**`annotate-stage`** (`0d0f517`, `eeeb6b6`; PR #31 `86c5677`, branch deleted 2026-11-06):
`--stage annotate`.

## Resolved decisions

- **`--stage` instead of `-entry` (2026-10-21):** Nextflow 26.04's strict parser rejects
  `-entry`; the entry workflow dispatches on a required `--stage`. Job scripts changed accordingly.
- **SingleR dropped (2026-10-21):** focus on CellTypist and scANVI instead.
- **`--qc_only` and `--prep_cellxgene` dropped (2026-10-22):** both were broken and unused.
- **Final h5ads always sparse + gzip (2026-10-22):** enforced in `ADATA_PUBLISH`; intermediates
  stay uncompressed.
- **Automatic filtering floors (2026-10-22):** `--min_genes` / `--min_counts` floor the N-MAD
  lower bounds; empty cells after gene filtering are always removed. Docs recommend low floors
  with automatic filtering.
- **D2, per-sample QC thresholds (2026-10-22):** implement; samplesheet value wins, empty cell
  falls back to the global parameter; keep the column name `min_counts_cell`; add `max_genes` /
  `max_counts`.
- **Unused parameters removed (2026-10-22):** `--mito_genes`, `--skip_enrichment` (both to the
  backlog for re-assessment) and `--symbol_col`.
- **D3, `exclude_samples_*` (2026-10-22):** implement, not remove (removing would not even make a
  differentialabundance contrasts file fail, since nf-schema ignores unknown columns).
- **`defaultBranch` (2026-09-25):** `dev_ukdri`, also GitHub's default branch (`bdad8e6`).
- **Container build route (2026-09-24):** the two remaining images were rebuilt from the committed
  Dockerfiles and published on Docker Hub: `docker.io/nhecker/scanpy-report:1.11.4-coreinf0.4`
  (report and enrichment) and `docker.io/nhecker/pydeseq2:0.1` (`DIFFERENTIAL_GENES_PER_CONTRAST`).
  Consequence for upstreaming: contributing `pydeseq2/differential_genes` to `nf-core/modules`
  means rebuilding it through Seqera Containers.
- **Container hosting (2026-09-01):** public Docker Hub is allowed for this repo's modules. Seqera
  Containers stays preferred where a plain `environment.yml` suffices, and is **required** for any
  module earmarked for `nf-core/modules`. GHCR was considered and rejected for the same reason.
  The choice is per module, not repo-wide.
