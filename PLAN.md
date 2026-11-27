# UK DRI scdownstream: clean up the fork, then port it to a new repo

> **Status 2026-11-23** (checked against `dev_ukdri` `bb5dd9f`, PR #34). All three stages
> (`qc_clustering` → `downstream` → `differential_genes`) complete on the UK DRI cluster with
> public images, reports included (SCR_MM_1 and SCR_MM_2 re-run on Nextflow 26.04). Since
> 2026-10-01: `--ambient_correction none` by default (PR #17), a differential_genes marker JSON
> (PR #18), Nextflow 26.04 strict-parser support with `--stage` replacing `-entry` (PR #19),
> SingleR removed (PR #20), sparse + gzip final h5ads (PR #21), `--prep_cellxgene` / `--qc_only`
> removed and `--duplicate_var_resolution max` fixed (PR #22), no empty cells after filtering and
> `--min_genes` / `--min_counts` as floors for automatic thresholds (PR #23), and per-contrast
> sample exclusion via `exclude_samples_col` / `exclude_samples_values` (PR #24, D3), the unused
> `--mito_genes` / `--skip_enrichment` removed (PR #25), and per-sample QC thresholds from the
> samplesheet with `--symbol_col` removed (PR #26, D2), and a clearer differential_genes report
> (PR #27, findings #44-#47). PR #28 fixed a LIANA+ dotplot crash in the downstream report
> found on SCR_HS_1, and PR #29 limited the reports' CellTypist labels to the top
> `--report_max_celltype_labels` (default 50). PR #30 (2026-11-06) added scANVI reference mapping in
> `qc_clustering` and `--stage train_models`, and PR #31 `--stage annotate`. PR #32 (2026-11-13) added presto/COSG marker genes and moved every main h5ad step to anndata 0.12.2 with own images. PR #33 (2026-11-23) left `--max_genes` / `--max_counts` unset by default, and PR #34 (2026-11-23) made Seurat v3 HVGs the default. What is done, and how, is in
> [PLAN_DONE.md](PLAN_DONE.md). This file lists only the open work. The open
> findings of the model-training feature are in
> [Feature: custom model training](#feature-custom-model-training-and-reference-mapping-celltypist-scanvi).

> **Current work (2026-11-23).** PR #33 (`--max_genes` / `--max_counts` unset by default, merge
> `8e32a86`) and PR #34 (Seurat v3 HVGs by default, merge `bb5dd9f`) are merged; both tested on
> the HPC first (see [PLAN_DONE.md](PLAN_DONE.md)). Before that, PR #32 (`anndata-0.12-throughout`, merge `20bb8f0`) was merged into
> `dev_ukdri`. It includes `marker-genes-presto-cosg` (presto/COSG markers), one anndata generation
> (0.12.2) for every main h5ad step, the own images `nhecker/scanpy-ukdri-suite:1.11.4-1` and
> `nhecker/seurat-ukdri-suite:5.5.1-1`, `ADATA_TORDS` / `.rds` outputs / SCimilarity / pseudobulking
> removed, `utils` subworkflow rename, and the downstream report legend and sortable tables; see
> [PLAN_DONE.md](PLAN_DONE.md). All five stages passed on real data (user, 2026-11-13).
>
> Open items:
> - HVG loose ends from PR #34 (not changed there): the samplesheet column `n_hvgs`
>   (assets/schema_input.json, docs/usage.md) has no effect since the dead `ext.n_hvgs` is gone;
>   QC_CLUSTER and INTEGRATE HVG outputs publish to the same folder with `--save_intermediates`
>   (pre-existing); the training report's HVG sentence also shows when only CellTypist is trained.
> - **Real-data run of the `seurat-ukdri-suite` modules** (not covered by the tests above):
>   `--ambient_correction decontx|soupx`, `--doublet_detection scds`, `--integration_methods seurat`.
>   The R templates (fixed for anndataR 1.2.1) only ran locally on PBMC data.
> - Regenerate nf-test snapshots: they still list `ADATA_TORDS`, `SCANPY_RANKGENESGROUPS`, and old
>   module versions (scvi/scanvi/celltypist/hvgs/add_obs_obsm, quality_control, filter, pipeline
>   tests). The nf-test CI fails until then. Needs a machine with Docker.
> - Step 3 of the anndata alignment: the nf-core modules in the per-sample QC path still have their
>   own anndata: `anndata/getsize` and `anndata/barcodes` (0.10.9, reads backed), `scanpy/scrublet`,
>   `doubletdetection`, `scvitools/solo|scar` (scvi-tools 1.2.0), `cellbender`. Nothing in the
>   pipeline writes `None` before them, but a user input h5ad with `null` elements fails there.
>   Patch them to the base/suite image. scvi-tools (scVI/scANVI) stays on Wave (anndata 0.12.1).
> - `SCANPY_ENRICH` and `SCANPY_EXPORT_MARKERS` have no `stub:` block; `ADATA_EXTEND` neither.
> - New normalisation methods (Pearson residuals, SCTransform with a merge subworkflow), see
>   Feature: normalisation methods.
> - Marker genes at large scale: memory first (gene chunks, no full CSC copy), then optional
>   parallelism; see Backlog: marker genes at large scale.
> - Container generation upgrade (pandas 3 / anndata 0.13), all images together; see Backlog:
>   move to pandas 3 / anndata 0.13.
> - Other "scdownstream" names (report templates `scdownstream_*.qmd`, unused `.ipynb` report
>   templates, logos, `--name` default): renamed at the port to `sc-analysis-nf` (Stage B, D1).
> - #53 (command-line params are Strings under Nextflow 26): typed `params {}` block in `main.nf`,
>   own branch; docs already recommend `-params-file`.
> - SYNCED_COMMIT: `bb5dd9f` (HPC synced 2026-11-23).

## Aim

The UK DRI fork of `nf-core/scdownstream` restructured a single-pass pipeline into three
sequential stages, selected with `--stage` (Nextflow 26.04's strict parser no longer supports
`-entry`). It now runs, but it still carries silent-correctness bugs, dead
parameters and dead code, in-place edits to nf-core modules, and almost no tests. The plan has
two stages:

1. **Stage A: a clean, working version in this fork.** All of the fixing and cleanup happens
   here: correctness fixes, dead code and dead parameters removed, per-stage parameter schemas,
   completion/versions, tests and CI. The fork becomes the reference implementation. It ends with
   the `v0.9.0` tag.
2. **Stage B: port to a new repo.** A mechanical copy, not a rewrite:
   - copy the non-nf-core modules (`modules/local/*`) as they are;
   - reinstall the nf-core modules (`modules/nf-core/*`) and register patches so they match the
     fork;
   - copy the workflows, subworkflows and entry points.

   Because the new repo only ports the fork, a run of it must give the same results as `v0.9.0`.
   Any difference is a porting bug.

This replaces the earlier plan, which re-authored and triaged modules in the new repo and kept the
fork frozen apart from blockers.

## Open decisions

| # | Decision | Blocks |
|---|---|---|
| D1 | ~~Pipeline name and GitHub org path for the new repo.~~ **Decided 2026-11-18: `UKDRI/sc-analysis-nf`** (`nextflow run ukdri/sc-analysis-nf`). Empty repo created by the user 2026-11-18 to reserve the name; nothing was to move until the HVG changes were tested (done 2026-11-23, PR #34). First step: `DESIGN.md` in the new repo. The `-nf` suffix marks a Nextflow pipeline (as nextflow-io's `rnaseq-nf`) without the `nf-` plugin prefix or an nf-core look-alike; no `ukdri` (the org carries it), no `scdownstream`/`nf-core`. Do **not** use the old placeholder `scflow`: nf-core/scflow was the UK DRI/Imperial pipeline, archived 2024. | — |
| D2 | ~~Per-sample QC columns (#2).~~ **Resolved 2026-10-22:** wired in; a samplesheet value overrides the global parameter, an empty cell falls back to it; `max_genes` / `max_counts` columns added (PR #26). | — |
| D3 | ~~`exclude_samples_*` (#17).~~ **Resolved 2026-10-22:** implemented per contrast in `deseq2_diffgenes_per_contrast.py`; the column must be sample-level, which pseudobulk keeps (PR #24). | — |
| D4 | ~~`--ambient_correction` default (#24).~~ **Resolved 2026-10-06:** default `none` (branch `ambient-default-none`). | — |
| D5 | ~~Scrublet unscorable cells get `NaN` → `True`.~~ **Resolved 2026-10-01:** they count as not doublets (`fillna(False)`, PR #16). | — |

## Open findings

Numbers are the original audit numbers, kept so commits and `PLAN_DONE.md` can refer to them.
Line references re-checked on 2026-09-29.

### Silent wrong results

**3. Clustering silently depends on the scVI neighbour graph.**
`params.cluster_neighbors` defaults to `neighbors_scvi` ([nextflow.config:66](nextflow.config#L66)),
produced only by [integrate.nf:96](subworkflows/local/integrate.nf#L96).
[qc_clustering.nf:162](workflows/qc_clustering.nf#L162) builds `neighbors_pca`, then
[line 170](workflows/qc_clustering.nf#L170) asks Leiden for `neighbors_scvi`. Dropping `scvi`
from `--integration_methods` gives no clustering and no error. (README limitation 2.)

**3a. A bare `except:` turns any failure into empty enrichment results.**
[enrich.py:63](modules/local/scanpy/enrich/templates/enrich.py#L63) wraps `sc.queries.enrich()`
in `try: … except:`, prints a warning and writes `pd.DataFrame({'query': ['empty set']})` for the
group. With `gprofiler-official` now in the image the missing-dependency case is gone, but a
network failure still gives an empty table for every group, a zero exit code and a report that
renders with nothing in it. Narrow the clause to the exceptions g:Profiler actually raises.

### Cannot run, or needs a workaround

**4c. Three conda envs are not reproducible.** `pydeseq2/differential_genes` and both
`decoupler/*` envs still list every dependency under `pip:` with unpinned `>=` bounds, although
`pydeseq2`, `decoupler`, `scanpy` and `anndata` are all on conda-forge/bioconda. No env has been
solved (`micromamba create --dry-run`).

**14. Two parameters are defined twice in `nextflow.config`**: `outdir` (17, 139) and `email`
(19, 141). Both copies are `null`, so nothing breaks yet. (`celltypist_model` was de-duplicated and
`celldex_reference` removed in PR #20.)

### Dead parameters and dead code

**16. Documented-as-non-functional parameters still appear in `--help`**: `--skip_liana`,
`--skip_rankgenesgroups`, `--pseudobulk*`, `--cluster_per_label`,
`--cluster_global`, `--unify_gene_symbols`, `--base_embeddings`, `--base_label_col`
(`--duplicate_var_resolution` and `--aggregate_isoforms` are live: both are used by `ADATA_UNIFY`)
(`--symbol_col` was removed in PR #26; gene symbols are per-sample `symbol_col` only)
(`--prep_cellxgene` and `--qc_only` were removed in PR #22; `--doublet_detection_threshold` is live since `--doublet_removal`, 2026-10-01). One flat parameter namespace across three entry points is why this
class keeps growing: a parameter meaningful to one stage is silently inert in another.

**17a. `expected_cells` in `assets/schema_input.json` has `minimum: 1` but `default: 0`**
(found 2026-10-22), so the default itself would fail validation if it were ever checked.

**18. Dead code**:
- `DIFFERENTIAL_GENES_PER_CLUSTER` (never included anywhere);
- `QC_FILTERED`, imported but commented out in `quality_control` (see #22). (`DOUBLET_REMOVAL`
  and `GET_DEDOUBLETED_SIZE` are in use since `--doublet_removal`, PR #16.)
- `unify_genes.nf` and `hugounifier/{get,apply}`, reachable only via the retired
  `--unify_gene_symbols`;
- the three `.ipynb` templates in `modules/local/scanpy/report/templates/`;
- five unused imports in [downstream.nf:8-12](workflows/downstream.nf#L8-L12) (`MULTIQC`,
  `paramsSummaryMap`, `paramsSummaryMultiqc`, `softwareVersionsToYAML`,
  `methodsDescriptionText`);
- dead `ch_label_grouping` / `grouping_col` in [qc_clustering.nf:150-151](workflows/qc_clustering.nf#L150-L151);
- since the default workflow was removed (`3aa34e3`): `subworkflows/local/pseudobulking`, which
  only it used (`cluster` and `per_group` were deleted in PR #19, the unused `FINALIZE` subworkflow
  in PR #22);
- the `tests/*.nf.test` pipeline tests, which target the removed default entry, so their
  snapshots are stale (`main_pipeline_qc` was deleted with `--qc_only` in PR #22).

**23. Intermediate results are passed between processes as pandas `.pkl` files** *(low
priority)*. 20 modules emit a `.pkl` (obs/var/obsm fragments from scrublet, solo,
doubletdetection, celltypist, hvgs, pca, umap, leiden, paga, rankgenesgroups, harmony, combat,
scvi/scanvi, scimilarity, liana, entropy, merge, mergeembeddings), and `adata/extend` unpickles
them. Pickle ties the writer's and reader's pandas/numpy versions together across containers,
can't be inspected without Python, and unpickling runs arbitrary code. See Backlog.

### Tests, CI and completeness

**9. Zero coverage for anything the fork added.** No test exists for `qc_clustering`,
`downstream`, `differential_genes`, or any of the eight added modules. Only `qc_clustering` has
test data (the upstream samplesheet). After PR #16 the snapshots of the `DOUBLET_REMOVAL`,
`DOUBLET_DETECTION` and `QUALITY_CONTROL` nf-tests are stale (new MultiQC text, method names, and
versions now listed for every doublet method) and need `--update-snapshot`.
The `SCANPY_FILTER` module tests were also out of date (8 inputs passed to an 11-input process);
their inputs were fixed when `--mito_genes` was removed (2026-10-22), but the snapshots need
regenerating, the "symbol column that does not exist" test expects a failure that can no longer
happen (`symbol_col` is not read by `filter.py`), and it is unchecked whether the test h5ad has the
QC columns (`total_counts`, `n_genes_by_counts`, `pct_counts_mt`) that `filter.py` reads.

**10. `.github/workflows/nf-test.yml` triggers only on `dev`/`master`/`main`**, so `dev_ukdri`
and the feature branches never run tests.

**12. `modules/nf-core/scanpy/scrublet` was edited in place with no patch registered** in
`modules.json`, unlike `cellbender/merge`, `scvitools/scar` and `scvitools/solo`, which carry
`.diff` files. `nf-core modules update` would silently clobber it, and Stage B B3 cannot reproduce
it. The edits are: the h5ad written after the `bool` cast (PR #7) and, since PR #16, unscorable
cells (`NaN`) counted as not doublets (`fillna(False)`).

**13. Identity is still upstream's** (new-repo work): `manifest.name = 'nf-core/scdownstream'`,
`.nf-core.yml` `is_nfcore: true`, `assets/schema_contrasts.json` `$id` → `nf-core/scdownstream/master`,
and the nf-core-only CI workflows (`branch.yml`, `awstest.yml`, `awsfulltest.yml`,
`release-announcements.yml`, `template-version-comment.yml`, `linting_comment.yml`).

**19. `downstream` and `differential_genes` never call `PIPELINE_COMPLETION`**
([main.nf:61](main.nf#L61) is the only call, in `qc_clustering`), so two of three stages have no
completion summary, no `--email` / `--hook_url` and no failure notification.

**20. Version capture is incomplete**: `qc_clustering` mixes versions from ~4 of ~15 processes;
`downstream` mixes only `FINALIZE_H5AD`; `downstream` and `differential_genes` never write a
`software_versions.yml`. (The report modules' `versions.yml` output was restored in `d28a17e`.)

**22. No MultiQC for `downstream` / `differential_genes`**, and no post-filtering QC plots
anywhere (`QC_FILTERED` is commented out).

*Found 2026-09-29 while stub-testing `--metadata` (branch `add_metadata`):*

**53. Command-line parameters are Strings under Nextflow 26's strict parser** (found 2026-10-28).
`--flag false` arrives as the String `"false"`, which Groovy treats as true; numbers arrive as
Strings too (`--n 200` → `"200"`). Parameters from `-params-file` YAML keep their types, so cluster
jobs built with params.yml are not affected. Affected: every Groovy-level use of a boolean or
numeric parameter given on the command line, e.g. a `--x false` overriding a `true` default, the
`--neighbors_n_pcs` / `--pca_n_comps` comparison in `validatePcaParameters()`,
`integrate.nf`'s `n_hvgs >= 0`, and `publishDir enabled:` / ext closures in the config.
Verified fix: a typed `params { name: Type = default }` block in `main.nf` (Nextflow 25.10+) makes
Nextflow convert command-line values (Boolean, Integer, Float). Migrate all parameters in a
separate branch; `train_models` already converts its own booleans with `paramAsBoolean()`.

**25. `--save_intermediates` crashes `qc_clustering`.** The `withName: 'GET_.*_SIZE'` block in
`conf/modules.config` sets `publishDir = [enabled: params.save_intermediates]` with no `path`, so
with the flag on `GET_FILTERED_SIZE` fails: "Target path for directive publishDir cannot be null".

**26. A pipeline `-stub` run cannot get through `qc_clustering`.** 41 local module stubs
`touch` an empty `versions.yml` (`ADATA_PUBLISH` was fixed in PR #21, `PYDESEQ2_GENERATE_REPORT` in PR #27), and `processVersionsFromYAML`
(`subworkflows/nf-core/utils_nfcore_pipeline/main.nf:82`) calls `collectEntries()` on the null
it parses, aborting at the first finished task. `ADATA_EXTEND` (`FINALIZE_QC_ANNDATAS`) has no
stub block, so it runs real Python in stub mode. Stubs should write a real `versions.yml`.

---

## Stage A: a clean, working fork

Each step is a small branch off `dev_ukdri`, reviewed and merged by PR, as so far.

### A1. Correctness

1. **#3:** make the Leiden neighbours key follow what the stage actually built, or stop with a
   clear error when `scvi` is missing from `--integration_methods`. Never produce nothing.
2. **#3a:** narrow the bare `except:` in `enrich.py`.
3. ~~Apply decisions D2 (#2) and D3 (#17).~~ Both done: D3 in PR #24, D2 in PR #26.

### A2. Cleanup

1. Remove the remaining duplicate params, `outdir` and `email` (#14).
2. Retire every dead parameter in #16. If a capability is wanted later, it comes back with an
   implementation and a test.
3. Delete the dead code in #18, including the stale `tests/*.nf.test`.
4. `QC_FILTERED`: re-enable it (post-filtering QC plots for MultiQC, #22) or delete the import.
   Not the commented-out middle ground. (`GET_DEDOUBLETED_SIZE` runs with `--doublet_removal`
   since PR #16.)
5. Give the `GET_.*_SIZE` publishDir a path, or disable it (#25).
6. Make every stub write a real `versions.yml` and give `ADATA_EXTEND` a stub (#26), so that
   `-profile test -stub` runs end to end. That is the cheapest local wiring test for every entry.

### A3. Per-stage parameter schemas

1. Keep the three entry points in one repo; they share the AnnData contract and most modules.
2. Give each entry point its own schema: `assets/schema_params_qc_clustering.json`,
   `..._downstream.json`, `..._differential_genes.json`. nf-schema's `validateParameters()`
   accepts a `parameters_schema:` path, so each entry validates against only its own parameters
   and rejects the others by name instead of ignoring them.
3. `--input` is required only by `qc_clustering`; the other two take `--base_adata`.
4. Keep `nextflow_schema.json` valid for `nf-core pipelines lint` (e.g. the union, or the
   `qc_clustering` schema), and say in the docs which schema each stage uses.

### A4. nf-core module hygiene (so Stage B is mechanical)

1. **#12:** register the scrublet edit with `nf-core modules patch scanpy/scrublet`.
2. Run `nf-core modules lint` / `nf-core modules patch` over every `modules/nf-core/*` and confirm
   each difference from its recorded `git_sha` is either a registered `.diff` or nothing.
3. For the UK DRI-authored modules (`decoupler/{pseudobulk,filter_pseudobulk}`,
   `pydeseq2/{differential_genes,report}`, `scanpy/{enrich,report,export_markers,normalization}`,
   `adata/{publish,add_obs_obsm}`), bring them to the nf-core/modules standard now, while it is
   cheap: a `meta.yml`, an nf-test with a stub run in the `tests/main.nf.test` layout, and no
   `params.*` inside the module. This keeps the generic ones (decoupler pair, pydeseq2)
   contributable later.

### A5. Reproducibility

1. **#4c:** rewrite the `pydeseq2/differential_genes` and `decoupler/*` envs as pinned conda
   dependencies, matching the versions in the published images. Resolve every bound against the
   API the template calls, not the newest release (the `decoupler>=1.6` / `dc.pp` mismatch is the
   failure mode).
2. Solve every env once: `micromamba create --dry-run -f <env>`.
3. Pin `gcfntnu/scanpy:1.11.4` by digest
   (`sha256:a97059ba131a56b9c0e7b8ddaf41fd7d4d781d0ebb2bd04759d7afc7973c0e44`, pushed
   2025-08-26), and the `nhecker/*` images too.

### A6. Completeness

1. `PIPELINE_COMPLETION` on all three entry points (#19).
2. Mix `versions` from every process and write `software_versions.yml` for every stage (#20).
3. Extend MultiQC to `downstream` and `differential_genes`, or state in the docs that the Quarto
   reports are the intended view for those stages (#22).

### A7. Tests and CI

1. One nf-test per stage, each with its own `--stage`: `qc_clustering` from the samplesheet,
   `downstream` and `differential_genes` from a small committed fixture `.h5ad` plus a
   `contrasts.tsv`, so the stages test independently and without private data (#9).
2. One nf-test per UK DRI-authored module (A4.3), following the pattern in
   `modules/local/*/tests/main.nf.test`.
3. An nf-test for `QUALITY_CONTROL` with `ambient_correction = 'none'` (#24).
4. CI (#10): run nf-test on `dev_ukdri` and PRs into it, with `-profile test,docker` and
   `-profile test,singularity`.
5. **Container-resolution check** (script first, CI job after): grep every `container` directive
   under `modules/`, resolve each manifest against its registry, and fail on anything that is not
   HTTP 200. Tokens: Wave via
   `https://cerbero.seqera.io/auth/token?service=community.wave.seqera.io&scope=repository:library/<name>:pull`,
   Docker Hub via `https://auth.docker.io/token`, quay.io needs none. ~54 requests, a few seconds.
   Also fail on the literal `CONTAINER_REGISTRY_MISSING`, and flag any Docker Hub reference pinned
   by a mutable tag rather than a digest. Run it on a schedule as well as per PR, so a withdrawn
   image shows up as a failed nightly, not a failed user run. Consider committing it with the
   Appendix A smoke test under `tests/containers/`.

### A8. Off-cluster proof

1. An off-cluster x86_64 Linux run of all three stages (a CI runner or a VM with no `/nfsdata`).
2. `-profile test,docker` and `-profile test,singularity` runs of all three stages.
3. `-profile conda` smoke test. Quarto-from-conda rendering is the one step not yet proven.

### A9. Freeze

Tag `dev_ukdri` as `v0.9.0` once A1–A8 pass the Stage A acceptance checks, and say in
`README.md` that this tag is the supported version until the new repo ships. After that the fork
gets only fixes for blockers users find.

## Feature: custom model training and reference mapping (CellTypist, scANVI)

*Phase 1 (scANVI reference mapping in `qc_clustering`) and Phase 2 (`--stage train_models`) are merged (PR #30, 2026-11-06), together with `--stage annotate` (PR #31); decisions, phases and verification are in [PLAN_DONE.md](PLAN_DONE.md). Open below: the scVI / scANVI / scArches findings and the later items.* Users want to train their own reference models and use them
to annotate new data, starting with CellTypist and scANVI, more methods later.

### Findings: scVI / scANVI / scArches (open)

The scVI/scANVI/scArches code is inherited from upstream and has never been used here. A code
review on 2026-09-29 (key lines verified) found that **reference mapping does not work today**.

*Blocking:*

**28. Supplied models are never loaded.** `subworkflows/local/integrate.nf:79,106` tag them
`id: 'scvi_model'` / `'scanvi_model'`, but `scvi.py:32-34` and `scanvi.py:28,38` check for
`'scvi'` / `'scanvi'`. A new model is trained from scratch without a warning. scANVI then either
fails ("Not enough labels to run scANVI") or trains on the query's own labels, so there is no
label transfer. Upstream bug; no test covers the reference-model path.

**29. scANVI results never reach the output.** `label:scANVI` and the scANVI latent travel in
`ch_obs` / `ch_obsm`, which `workflows/qc_clustering.nf:145` collects and never uses.
`FINALIZE_H5AD_SIMPLE(ch_h5ad)` (`:183`) only gets the scVI-based object.

**30. No output without scVI.** `INTEGRATE.out.h5ad` is only set in the scVI branch
(`integrate.nf:93-99`), so `--integration_methods scanvi` alone silently produces nothing after
integration.

**31. `--base_adata` reference mapping is broken.** `merge_embeddings.py:15,23` needs
`X_scanvi` and `label:scANVI` in the base, which no output has (KeyError). The output holds query
cells only (`combine.nf:40,52-66`). Once #28 is fixed, scVI would concatenate the base latent with
a latent from a separately trained query model. The base+query object also mixes log-normalised
and raw `X` and loses `counts` (`merge.py:61`). Query cells whose names already exist in the base
are dropped silently (`merge.py:77`).

*Wrong or misleading results (remain after the blockers are fixed):*

**32. scANVI is fitted to log-normalised data.** `SCANVI.setup_anndata` has no
`layer="counts"` (`scanvi.py:51`), unlike scVI, although both use a count likelihood (ZINB). The
query path also lacks a `counts` layer, so a model built with `from_scvi_model` fails on it.

**33. Query labels are not hidden.** They are not set to the model's unlabelled category before
`load_query_data` (`scanvi.py:28-30`). Labels unknown to the model raise
`Category … not found in source registry` (scvi-tools 1.3.3 does not extend label categories for
query data); known labels bias the fine-tuning, so `predict()` partly echoes the input labels.

**34. Query fine-tuning drifts the reference and uses the wrong epochs.** The default
`weight_decay` (1e-6) still shrinks the frozen reference weights; scArches uses
`plan_kwargs={"weight_decay": 0.0}`. `--scvi_max_epochs` drives reference training,
`from_scvi_model` fine-tuning and query mapping at once; without it, a query falls back to up to
400 epochs. It needs a separate query setting (scArches: ~100).

**35. Hard labels only.** `predict()` without `soft=True` gives no confidence, so confident and
forced transfers (e.g. cell types missing from the reference) look the same.

**36. The query gets its own HVGs.** With only `--scanvi_model`, scANVI's query input is the
query's own HVG subset (`integrate.nf:104`), so `prepare_query_anndata` zero-pads most reference
genes; only a buried 0.8-overlap warning. The query should keep all genes (counts).

**37. `UNIFY` rewrites labels, and distinct labels can collide.** `to_Florent_case`
(`unify.py:85-108`, applied at `:184`) drops every character other than letters, digits and `_`,
removes one trailing "s", and `capitalize()` lower-cases everything after the first letter:
"OPCs" → `Opc`, "Excitatory neurons L2-3" → `Excitatory_neurons_l2_3`, and "CD4+ T cells" and
"CD4- T cells" both → `Cd4_t_cell`. "nan"/"na"/"null"/"unknown" correctly become `unknown`.
Also at `unify.py:176-179`: `"unknown" in adata.obs["label"]` checks the index, not the values,
and `adata.obs["label"].replace({...}, inplace=True)` on a column taken from `obs` likely does
nothing under pandas copy-on-write, so a custom `unknown_label` would stay a class (to verify).
Scope today: only the samplesheet `label_col` in `qc_clustering` (not CellTypist/SingleR
predictions or `--metadata` columns). Not blocking for training, which gets its own prep step;
a `qc_clustering` fix for later.

*Not ideal:*

**38.** Several smaller issues:
- The model must be a single `model.pt` (`stageAs: 'reference_model/model.pt'`); a model
  directory fails. It is passed as a string into a `path` input rather than via `file(...)`.
- The scVI model is not published without `--save_intermediates`.
- Only `X_emb` is written (no `X_scanvi`), and there are no scANVI neighbours or UMAP.
- HVGs use `seurat` on log data without a batch key; a reference model needs `seurat_v3` on
  counts, per batch.
- Values are hard-coded in calls: `"label"`, `"unknown"`, `"batch"`, the seed,
  `early_stopping=True`, and matmul precision `"medium"` (bf16 on GPU). A comment says "Round to
  ensure hashes are stable", but no rounding follows.
- The parameter checks miss `--scanvi_model` without `scanvi` in `--integration_methods`, and
  `scanvi` without `scvi`.
- The "Reference mapping" docs (`docs/usage.md:559-574`) describe behaviour that does not
  happen.

**39. Label consumers assume CellTypist.**
- Both reports pick label columns with `startswith('celltypist:')`: the QC report at `:317`, the
  analysis report at `:92`, and its overlap table at `:147-210`.
- The analysis report casts the `--selected_clustering` groups to `int` (`:138`), so a cell type
  column as grouping fails. `--selected_clustering` and `--diffgenes_group_col` themselves
  accept any column.
- Column naming is inconsistent across tools: `celltypist:<model>` (+`:conf`), `label:scANVI`
  (no confidence), and `annotation:scimilarity` (SingleR, whose names included the sample id, was removed in PR #20).
- `CELLTYPES_CELLTYPIST` only accepts built-in model names: the schema pattern `[a-zA-Z0-9_]`
  rejects paths, and the template always calls `download_models` (`celltypist.py:54-58`).

### Later

- **HVG flavour for training (decided 2026-10-28):** `--stage train_models` selects the scANVI
  reference HVGs with scanpy's count-based `pearson_residuals` (`sc.experimental.pp.
  highly_variable_genes`, batch-aware). The planned `seurat_v3` needs `scikit-misc`, which the
  `SCANPY_HVGS` container (Wave `pyyaml_scanpy`) lacks. Re-evaluate whether `seurat_v3` is the
  better choice; if so, build a custom container with scanpy + scikit-misc and switch
  `ext.flavor` for `HVGS_EVAL` / `HVGS_FINAL`.
- `--base_adata` reference mapping (#31): emit the base+query object, align `X`/counts, and
  require a base produced by `train_models`.
- scANVI neighbours / UMAP.
- Consistent SingleR / SCimilarity column names.
- A reference latent UMAP in the training report.
- More training methods.

## Feature: normalisation methods (Pearson residuals, SCTransform)

*Added 2026-11-10, not started. Asked for by the user; nothing in this plan covered it before (only
Pearson residuals as the HVG flavour of `--stage train_models`, see Feature → Later).*

### Current state
- One method: `SCANPY_LOG_NORMALIZE` (`modules/local/scanpy/normalization`) runs
  `sc.pp.normalize_total` (median depth) + `sc.pp.log1p`; raw counts stay in `layers["counts"]`.
  Called in `qc_clustering.nf` and `subworkflows/local/integrate.nf`.
- The log-normalised `X` feeds HVGs, PCA, neighbours/UMAP/Leiden, marker genes (presto, COSG) and
  the report plots. scVI/scANVI, CellTypist, `--stage annotate` and pseudobulk read
  `layers["counts"]` and are not affected by the choice.
- SCTransform runs today only inside `SEURAT_INTEGRATION` (per batch, as part of the Seurat
  integration). That module converts the AnnData to Seurat with anndataR (`as_Seurat`) and writes
  back only `obsm["X_emb"]`.

### Goal
`--normalization log1p | pearson_residuals | sctransform` (default `log1p`, unchanged results).
- **Pearson residuals** (analytic, Lause et al. 2021): `sc.experimental.pp.normalize_pearson_residuals`
  (`theta`, `clip`), together with the Pearson-residual HVGs
  (`sc.experimental.pp.highly_variable_genes(flavor="pearson_residuals")`). The residuals are dense:
  compute them on the HVGs only and use them as the PCA input. They are not suitable for marker genes
  or expression plots, so `X` stays log-normalised for those.
- **SCTransform** (Seurat, `vst.flavor = "v2"` with glmGamPoi; Hafemeister & Satija 2019, Choudhary &
  Satija 2022): R module. Seurat returns corrected counts (`SCT` assay `counts`), their log1p
  (`data`), residuals for the variable features (`scale.data`) and the variable features.

### Proposed design
- **SCTransform as a subworkflow** (`NORMALIZE_SCTRANSFORM`), because anndataR may not keep
  everything when converting to Seurat and back (obsm, obsp, uns, varm, extra layers, categorical
  dtypes):
  1. `SEURAT_SCTRANSFORM` (R; Seurat + glmGamPoi + anndataR container) reads the h5ad, runs
     SCTransform (per sample/batch when given), and writes ONLY the new parts to a small h5ad:
     corrected counts, `data`, the residuals of the variable features, a `var` flag for the variable
     features, and the model parameters.
  2. A Python merge step (extend `ADATA_EXTEND` / `ADATA_ADD_OBS_OBSM`, or a new
     `ADATA_ADD_LAYERS`) adds these to the unchanged original object by cell and gene name, and checks
     that cells and genes match: layers `sct_counts` / `sct_data`, the residuals (variable features
     only, e.g. `obsm` or a separate HVG object for the PCA), `var["sct_variable"]`.
  `layers["counts"]` stays the raw counts.
- **Pearson residuals**: Python only (scanpy), in the normalisation / HVG / PCA steps; no
  conversion needed.
- **Which matrix where** (to decide, see below): the PCA input (residuals for both new methods), the
  HVGs (SCT variable features vs scanpy's), and `X` for markers and plots (`log1p` of
  `normalize_total`, or SCT `data` for `sctransform`).

### Open decisions
- `X` with `sctransform`: SCT `data` (Seurat's recommendation for visualisation; markers would then
  need `PrepSCTFindMarkers`-like handling) or the usual log-normalised counts.
- Per-sample/batch SCTransform (as in `SEURAT_INTEGRATION`) vs one model for all cells.
- HVG selection with `sctransform`: SCT variable features, or keep scanpy's HVG step.
- Container: check whether `docker.io/nicotru/seurat` has glmGamPoi and a current anndataR, or build
  one (pinned Seurat, glmGamPoi, anndataR).
- Memory: dense residuals and SCTransform on several hundred thousand cells (time, RAM).

### Verification
- SCTransform: compare the module's outputs with Seurat run directly in R on the same cells.
- Merge: every original slot (obs, var, obsm, obsp, uns, layers) unchanged after the merge; new
  layers aligned by name.
- Pearson residuals: equal to `sc.experimental.pp.normalize_pearson_residuals` on the HVGs.
- `qc_clustering` stub runs and a small real run with each method; `downstream` (markers) unchanged
  because `X` stays log-normalised (unless decided otherwise).

## Backlog: marker genes at large scale (memory first)

*Recorded 2026-11-10 from benchmarks on `marker-genes-presto-cosg`; not started. Priority set by the
user: memory before wall time.*

### Measurements (50k cells x 20k genes, 70M non-zeros, 20 groups, one Mac)
- **R presto** `wilcoxauc`: 3.5 s, 3.6 GB peak (including loading the data in R).
- **Our `MARKERS_PRESTO` template:** 17.7 s in total, of which about 16.5 s is the computation; 1.8 GB
  peak.
  - **Single-threaded:** 14.2 s user for 17.8 s wall. The per-gene loop runs in Python
    (`rankdata`, `bincount`), and `threadpool_limits` doesn't help there.
  - **Peak about 3x the matrix size:** the full CSC copy plus the boolean copy inside
    `sc.get.aggregate(count_nonzero)`.
- **Extrapolation to 5M cells** (about 2,500 detected genes per cell, about 12.5 billion non-zeros):
  - Python about 45-90 min, R about 15-25 min.
  - Memory decides: the matrix alone is 100-150 GB (64-bit indices), the Python peak 300-450 GB.
  - R's `dgCMatrix` can't hold more than 2^31 non-zeros at all.

### Step 1: memory (identical results, verified against R presto)
- **Gene chunks:** process genes in chunks (`X[:, chunk].tocsc()`, e.g. 1,000 genes) instead of a
  full CSC copy. That costs about 20 extra passes over the data, a few seconds.
- **Group statistics in the loop:** compute the group sums and non-zero counts inside the chunk loop
  instead of with `sc.get.aggregate`, to avoid its boolean copy. This is a documented exception to
  the "use library functions" agent rule; the reason is memory.
- **Target:** peak about 1.2x the matrix size; at 5M cells roughly 130-180 GB.

### Step 2 (optional): parallelism on top of the chunks
- numba `prange` / nogil threads, which share the matrix and also remove the per-gene Python
  overhead; or
- joblib with care: the `loky` backend memmaps large arrays into `/dev/shm`, which is RAM, i.e. one
  more copy of the matrix. Use a temp folder on disk or the fork backend, and expect about 150 MB per
  worker process.
- Expected roughly linear speed-up: the benchmark at about 1 s with 32 cores; 5M cells at about
  2-4 min.

### Step 3 (only if one node can't hold the matrix)
- Read gene chunks directly from the h5ad (backed mode). `X` is CSR on disk, so column reads are
  expensive; consider rewriting it once as CSC, or reading row blocks and accumulating.
- Check COSG's memory and core use the same way. Its authors report about 2 min for 1M cells, but
  it also needs the whole matrix.

### Not chosen
- R presto as the implementation: it needs an R container, the h5ad <-> R conversion costs time and
  memory and can lose parts of the object, and the `dgCMatrix` 2^31 non-zero limit rules it out at
  this scale.

## Backlog: move to pandas 3 / anndata 0.13 (all images together)

*Recorded 2026-11-10 after the COSG / SCANPY_ENRICH failure; not started.*

### Finding
- **The failure:** `MARKERS_COSG` ran in the BioContainer `cosg:1.2.1`, built with unpinned
  dependencies: pandas 3.1.1, anndata 0.13.4, scanpy 1.12.5. It wrote `obs_names` as a nullable
  string array. `SCANPY_ENRICH` (`nhecker/scanpy-report`, anndata 0.12.2) then failed on write:
  `allow_write_nullable_strings is False ... /obs`.
- **Reproduced locally.** With the pipeline's versions, the unchanged template works.
- **Fixed in `c973dc1`:** an own image, `nhecker/scanpy-cosg:1.11.4-cosg1.2.1`
  (`FROM gcfntnu/scanpy:1.11.4` + `pip install --no-deps cosg==1.2.1`).
- **The rule:** "Development gotchas" in `CLAUDE.md`, the agent rule, and the note in
  `docs/usage.md` "Container images".
- **Direction of failure:** mixing fails when a step with the older anndata has to *write* an h5ad
  that a new-generation step wrote. Reading works. Older-generation output read by a new step is
  fine.
- **A newer gcfntnu base doesn't help:** `gcfntnu/scanpy:1.12.3` (2026-09-03) was built with
  `pip install scanpy==1.12.3` and no further pins. scanpy 1.12.3 only requires `pandas>=2.3` and
  `anndata>=0.10.8`, and on that date pip resolved pandas 3.0.x / anndata 0.13.x. So it is the new
  generation too: switching single modules to it gives the same failure.

### State 2026-11-11 (after `anndata-0.12-throughout`)
- **Releases:** anndata 0.13.4 (2026-09-18; the 0.12 line ended with 0.12.19, 2026-07-02), scanpy
  1.12.5 (2026-10-29), pandas 3.1.1 (2026-11-04). anndata 0.13 and scanpy 1.12 need python >= 3.12,
  pandas >= 2.3, and anndata 0.13 needs zarr >= 3.1. The pipeline is on anndata 0.12.2 / scanpy
  1.11.4 / pandas 2.3.2 / python 3.11 (`gcfntnu/scanpy:1.11.4`, built 2025-08-26).
- **Much smaller job now:** after `anndata-0.12-throughout` the images are the base, our two suites
  (`scanpy-ukdri-suite`, `seurat-ukdri-suite`), scvi-tools (Wave) and the nf-core QC modules. Plan:
  an own pinned base (python 3.12, scanpy 1.12.x, anndata 0.13.x, pandas 3.x) instead of the
  unpinned gcfntnu builds, rebuild both suites on it, then scvi-tools >= 1.4 and the nf-core modules.
- **GPU images (checked 2026-11-11):** both suites are CPU only (no torch/CUDA). GPU steps use
  separate images: scVI/scANVI/solo switch under `-profile gpu` from Wave `scvi-tools:1.3.3` (CPU
  torch 2.7.1, python 3.13, anndata 0.12.1) to `ghcr.io/scverse/scvi-tools:py3.12-cu12-1.3.3-`
  (CUDA 12, built 2025-07-23 unpinned, so anndata 0.12.0/0.12.1 by date, not verified inside);
  scAR (scvi-tools 1.2.0) has no GPU variant; CellBender uses `broad-dsde-methods/cellbender:0.3.2`.
  With the migration: an own `scvi-ukdri-suite` on the new base with CUDA torch (runs on CPU too),
  one image for both profiles, so CPU and GPU runs use the same versions.
- **Side effects:** with anndata 0.13 writing nullable strings is allowed, so the anndataR pin
  (1.2.1, see `containers/README.md`) can be lifted; the R side must then read nullable strings
  (anndataR >= 1.2.2 does). Every tool must support python 3.12 (celltypist, liana < 3.14,
  pydeseq2, cosg, hugo-unifier: OK on PyPI metadata).

### Work (all at once)
- **Inventory:** list every image with its anndata/pandas/numpy versions. About 30 today:
  - Wave `pyyaml_scanpy` / `anndata_pyyaml` / `anndata:0.10.9`, `scvi-tools:1.3.3`,
    `celltypist_scanpy`, harmony / bbknn / leiden images;
  - the BioContainer `multiqc` (`liana:1.4.0` replaced by the report image, liana 1.6.1, 2026-11-11);
  - the Seurat/anndata2ri image;
  - `gcfntnu/scanpy:1.11.4`, `nhecker/scanpy-report`, `nhecker/pydeseq2`, `nhecker/scanpy-cosg`.
- **New or matching images** for every module that reads or writes h5ad. Our Dockerfiles move to
  the new base. Check the R side (anndataR / zellkonverter / anndata2ri) can read the new encodings.
- **Templates:** check for pandas 3 changes (default string dtype, copy-on-write, removed APIs) and
  the anndata 0.13 changes.
- **Pins:** update the `environment.yml` pins.
- **Tests:** end-to-end runs of every stage (qc_clustering, downstream, differential_genes,
  train_models, annotate) on the test data and one real dataset. Regenerate the nf-test snapshots.
- **Docs:** update the "one generation" note in `CLAUDE.md`, the agent rule and `docs/usage.md`.

## Feature: doublet detection clean-up

*Recorded 2026-10-01, not started.* All doublet methods run at one point, inside
`QUALITY_CONTROL` (QC metrics → `DOUBLET_DETECTION` → ambient correction → `SCANPY_FILTER`,
`subworkflows/local/quality_control/main.nf:68`), although they need different inputs. Only
scrublet is in the curated tool set today; the test profile also runs solo and scds.

### Findings (open)

**40. All methods run at one stage although they need different inputs.**
- **scrublet** (and scds, doubletdetection): scrublet simulates doublets by summing random pairs
  of the sample's own cells and scores observed and simulated cells by their PCA neighbourhood.
  These methods belong per sample on counts, after empty-droplet removal but **before** QC
  thresholding: upper thresholds (`max_counts`, `max_genes`) mostly remove real doublets and shift
  the score distribution. The current placement fits them. scrublet's `expected_doublet_rate` is
  left at its default of 0.06 (`scrublet.py:32`), although the rate scales with the number of
  cells loaded.
- **solo** is a classifier on an scVI model. scvi-tools recommends
  `SOLO.from_scvi_model(scvi_model, restrict_to_batch=<sample>)` on an scVI model trained on
  **QC-filtered** cells, typically the integration model, before clustering.
  `modules/nf-core/scvitools/solo/templates/solo.py:52-56` instead trains a new scVI model per
  sample on unfiltered cells and all genes, before filtering: slow, and the latent space also
  fits the low-quality cells.

**41. solo loop overwrites the scVI model.** In the per-batch loop
(`solo.py:64`), `model = SOLO.from_scvi_model(model, ...)` replaces the scVI model with the
SOLO model, so with more than one batch per sample the second iteration calls
`from_scvi_model` on a SOLO model.

**42. Only scrublet's call reaches `obs`** (`predicted_doublet`). solo, scds and
doubletdetection only write `.pkl`/`.csv` files, which only `--doublet_removal` reads. With
several methods, the annotation shows scrublet alone while removal uses the consensus.

**43. `--doublet_detection_threshold` mixes votes** from methods that belong at different stages.

### Proposed design (to decide)

Keep subworkflows, **split by stage**, rather than moving the methods into `qc_clustering.nf`.
Inlining would not fix the timing, would add one branch per method to the workflow, and would lose
the subworkflow nf-tests.

| Stage | Where | Methods |
|---|---|---|
| Per sample, before QC filtering | `DOUBLET_DETECTION_PREFILTER` in `QUALITY_CONTROL` (today's subworkflow, slimmed down) | scrublet, scds, doubletdetection; **scDblFinder (future)** |
| Per sample, after QC filtering | `DOUBLET_DETECTION_POSTFILTER` after `SCANPY_FILTER` | **DoubletFinder (future)**: needs a preprocessed, QC-filtered Seurat object per sample, a pK sweep and an expected doublet rate |
| After integration, before clustering | in `INTEGRATE` after `SCVITOOLS_SCVI`, or a small `DOUBLET_DETECTION_SCVI` | solo with `from_scvi_model(..., restrict_to_batch=sample)` on the integration model |

- **One list, routed by stage:** `--doublet_detection` stays a single list, and each method is
  sent to its stage.
- **Calls in `obs`:** every method writes `doublet:<method>` and `doublet:<method>:score`, and the
  QC report shows how the methods agree.
- **Doublet rate as a parameter:** `expected_doublet_rate` becomes a parameter.
- **Removal per stage:** each stage applies its own removal, instead of one mixed threshold.

### Open decisions

- Which methods stay in the curated tool set.
- How `--doublet_removal` / `--doublet_detection_threshold` work across stages. Options: a vote
  per stage, or annotate everywhere and remove only before filtering and at the solo step.
- Whether solo removal after integration is wanted, since removing cells after integration means
  the latent space was fitted with them included.

### Future methods

- **scDblFinder** (R/Bioconductor): pre-filter group, per sample on counts (`samples=` for
  multi-sample objects).
- **DoubletFinder** (R/Seurat): post-filter group.
- Both need containers, pinned `environment.yml` files and nf-tests, following the existing
  `scds` module (R template via anndata2ri).

*Decided 2026-10-01:* no multiplexing-based doublet detection (cell hashing or genotype
demultiplexing, e.g. HTODemux, vireo, souporcell). It needs HTO or genotype data that is rarely
produced together with scRNA-seq, and adds too much processing. The pipeline stays with
expression-based methods. When choosing the curated set, consider dropping scds and
doubletdetection rather than adding more expression-based methods.

## Feature: gene symbols from a GTF for inputs without symbols

*Recorded 2026-10-22, not started.* Some inputs, notably CellBender `.h5`/`.h5ad` output, lack the
var/obs annotation columns, including any gene-name column, so `var_names` are Ensembl IDs. Gene
symbols are needed for LIANA+ (ligand–receptor resources are symbol-based), enrichment, CellTypist
models and the mitochondrial QC metric.

### Findings (open)

**48. Without symbols, mito QC silently does nothing.** QC_RAW
(`modules/local/scanpy/plotqc/templates/plotqc.py:24-39`) looks for a `gene_name`, `gene_names`,
`gene_symbol` or `gene_symbols` column and otherwise flags mito genes with
`var_names.str.startswith("mt-")`. With Ensembl IDs that matches nothing, so `pct_counts_mt` is 0 for
every cell. `--max_mito_percentage` and the automatic mito threshold then remove nothing, and no
error is raised. QC_RAW also guesses symbol columns instead of using the samplesheet `symbol_col`.

**49. The existing fallback runs too late and needs the internet.** With `symbol_col = none`,
`UNIFY` (`subworkflows/local/unify.nf:16-25`) calls `ADATA_MYGENE`, an online lookup on
mygene.info. It runs after QC metrics, filtering and CellTypist, which all needed symbols already,
and compute nodes on the cluster may have no internet access.

**50. `symbol_col = none` crashes CellTypist, and MyGene is human-only** (found 2026-10-22).
CellTypist runs before UNIFY's MyGene conversion and treats `none` as a column name
(`celltypist.py:47-49`), so `symbol_col = none` with `--celltypist_model` stops with "Symbol column
none not found". `mygene.py:26` hardcodes `species="human"`. Both go away if symbols are resolved at
load time (this feature); until then, `docs/usage.md` ("Gene symbols") documents the limits.

### Proposed design

- **New parameter `--gtf`** (path, `.gtf` or `.gtf.gz`; schema `format: file-path`, `exists: true`).
  It should match the annotation used for quantification (e.g. the GENCODE/Ensembl release from
  nf-core/scrnaseq).
- **New module** (e.g. `ADATA_GTFSYMBOLS`):
  - parse the `gene` records of the GTF into a `gene_id` → `gene_name` map, with versions stripped
    on both sides (`ENSG00000141510.12` → `ENSG00000141510`);
  - keep the original IDs in `var["gene_ids"]`;
  - set the symbols as `var_names`, falling back to the Ensembl ID for genes without a `gene_name`;
  - **call `var_names_make_unique()`**: GTFs map several IDs to the same symbol (e.g. PAR genes,
    readthroughs), and duplicate `var_names` break several steps;
  - write the counts of mapped, unmapped and made-unique genes to the log and to `uns`.
- **Placement: before QC metrics, not in UNIFY.** QC_RAW's mito/ribo/hb flags and CellTypist need
  symbols (finding 48), so the mapping must run in `LOAD_H5AD` right after reading, before
  `ADATA_ADDSAMPLE` / QC_RAW. In UNIFY it would come after QC had already been computed on IDs.
  The only cost of running it early is one extra pass per sample.
- **Which samples:** apply the mapping only to samples whose `var_names` look like Ensembl IDs
  (`^ENS[A-Z]*G\d+`) and that have no usable symbol column. Samples that already have symbols are
  left untouched.

### Open decisions

- **`make_unique` vs `--duplicate_var_resolution`:** unifying early with `var_names_make_unique()`
  means UNIFY's `sum`/`mean`/`max` aggregation never sees these duplicates (they become `GENE-1`).
  The alternative is to keep duplicates until UNIFY and let `--duplicate_var_resolution` handle
  them, but QC_RAW and CellTypist would then run with duplicate `var_names`. The user asked for
  `make_unique`; confirm that the aggregation option should not apply to GTF-mapped duplicates.
- **Global or per sample:** a single `--gtf` for all samples, or also a samplesheet `gtf` column for
  mixed annotations.
- **ADATA_MYGENE:** retire it in favour of the GTF route (it works offline and comes earlier), or
  keep it as a fallback when no `--gtf` is given.
- **QC_RAW:** whether to also make it use the samplesheet `symbol_col` instead of guessing column
  names (finding 48).

### Verification

- Unit test of the module on a small GTF (versioned IDs, a missing `gene_name`, two IDs with the same
  symbol): check the map, `var["gene_ids"]`, unique `var_names` and the fallback.
- A run on a CellBender sample with only Ensembl IDs: `pct_counts_mt` is non-zero after the change,
  CellTypist runs on symbols, and LIANA+ finds ligand–receptor pairs.

## Stage B: port to the new repo

0. **`DESIGN.md` first** (decided 2026-11-18): before any code is moved, write a `DESIGN.md` in
   `UKDRI/sc-analysis-nf` that sets the rules for porting code from this fork (what is copied as is,
   what is rewritten or dropped, naming, containers, tests and docs required per module). Every
   later step follows it.
1. **Scaffold** with `nf-core pipelines create`, `is_nfcore: false`, org `ukdri`, skipping
   `fastqc` and `igenomes` as the current `.nf-core.yml` does. This lints against the non-nf-core
   ruleset, so UK DRI changes stop fighting template-sync. Keep only `nf-test.yml` and
   `linting.yml` from CI, triggered on `dev`, `main` and `feature/**`. Set `manifest` to the new
   name `ukdri/sc-analysis-nf` (D1), homepage, `defaultBranch = 'dev'` and a real semver (`1.0.0dev`), and re-point the
   `assets/schema_contrasts.json` `$id` (#13).
2. **Copy the non-nf-core modules**: `modules/local/*` as they are, with templates, envs,
   Dockerfiles, `meta.yml` and tests.
3. **nf-core modules**: `nf-core modules install` each module in the fork's `modules.json`
   (`anndata/{barcodes,getsize}`, `cellbender/{merge,removebackground}`, `doubletdetection`,
   `multiqc`, `scanpy/scrublet`, `scvitools/{scar,solo}`, `untar`), at the recorded `git_sha` or at
   the current upstream version, decided per module. Then re-apply the fork's changes and register
   them with `nf-core modules patch`: `cellbender/merge`, `scvitools/scar`, `scvitools/solo`,
   `scanpy/scrublet`. Unchanged modules stay unpatched. Also copy the nf-core subworkflows
   (`subworkflows/nf-core/*`) the same way.
4. **Copy the pipeline**: the `main.nf` entry points, `workflows/`, `subworkflows/local/`,
   `conf/`, the `assets/` schemas (samplesheet, contrasts, the three per-stage schemas), the test
   fixtures and nf-tests, and the container-resolution script. Rename the
   `utils` subworkflow to the new name. Rename the remaining "scdownstream" names: the report
   templates `scdownstream_*.qmd` (drop the unused `.ipynb` templates), logos and assets, the
   `--name` default, and the pipeline title in the reports and docs.
5. **Docs**: port `README.md`, `docs/usage.md` and `docs/output.md` rather than rewriting them.
   Add a migration note from `v0.9.0` (changed parameter names, and the changed
   `--ambient_correction` default if D4 goes that way). Keep `ACKNOWLEDGEMENTS.md` and
   `CITATIONS.md`; upstream attribution matters more in a separately named pipeline. Keep the
   docs free of development history: the reason behind each container choice
   belong here, in `PLAN_DONE.md` or in code comments.
6. **Later, optional**: contribute the generic modules (`decoupler/pseudobulk`,
   `decoupler/filter_pseudobulk`, `pydeseq2/differential_genes`) to `nf-core/modules`. That needs
   Seqera Containers images, not the `nhecker` Docker Hub ones (see Appendix B).

## Backlog (re-assess): custom mitochondrial gene list (`--mito_genes`)

*Removed 2026-10-22 (PR #25).* `--mito_genes` took a file of mitochondrial gene
names, was plumbed into `SCANPY_FILTER` as a `path` input, and was never used: the code reading it
was commented out in `filter.py`, and QC_RAW always flags mito genes by the `mt-` prefix. Removed
rather than implemented because the pipeline already assumes standard gene symbols elsewhere
(CellTypist, LIANA+, enrichment; see the GTF feature, #48), and with standard symbols the `mt-`
prefix is right for human and mouse. A custom list would only matter for unusual annotations.

Re-assess if a real dataset needs it (e.g. a species or annotation whose mitochondrial genes don't
carry the `mt-`/`MT-` prefix). If so, bring it back by computing `pct_counts_mt` in the filter
module from the list (falling back to the prefix), with a regression test whose retained cell
counts differ from the prefix default.

## Backlog (re-assess): an option to skip gene set enrichment (`--skip_enrichment`)

*Removed 2026-10-22 (PR #25), finding #15.* `skip_enrichment` was in
`nextflow.config` but not in the schema, and nothing read it: `SCANPY_ENRICH` always runs in
`downstream` (`workflows/downstream.nf`). Re-assess whether a working skip option is wanted.
The likely reason would be offline compute nodes: enrichment queries g:Profiler over the internet,
and a network failure currently gives empty tables without an error (#3a). If it comes back, it
needs a schema entry, a branch around `SCANPY_ENRICH` in `downstream.nf`, a report that copes with
missing enrichment results, and a test.

## Backlog (low priority): replace `.pkl` hand-offs with AnnData (#23)

Not a release blocker. Do it in the fork after Stage A, or in the new repo after Stage B.

1. Each module that emits an obs/var/obsm `.pkl` passes its results as AnnData instead: its
   `.h5ad` output directly, or a slim `.h5ad` holding only the new columns/embeddings (empty `X`).
2. Rewrite `adata/extend` (and the finalize subworkflows) to merge from those files, reading them
   with `anndata.read_h5ad(..., backed="r")` or only the needed element with
   `anndata.io.read_elem`, so the large base matrix is never loaded to copy a few columns.
3. Where a module only needs metadata from a large input, use backed mode too (the qc_clustering
   report already does, `075e9c9`).
4. Remove the `.pkl` outputs from `main.nf`/`meta.yml`, and update the nf-test snapshots.

Verify that the finalized `.h5ad` has identical obs/var/obsm to one built via the `.pkl` path.

---

## Verification

**Stage A acceptance (before the `v0.9.0` tag):**

1. **Every container reference resolves:** the A7.5 check reports zero misses and no
   `CONTAINER_REGISTRY_MISSING`. Run it first; it is the fastest signal.
2. **Portability:** on a machine with no access to `/nfsdata`,
   `nextflow run . -profile test,docker --stage qc_clustering --outdir out/`, and the same for
   `downstream` and `differential_genes`, all complete with no site path referenced. Repeat with
   `-profile test,singularity`.
3. **Parameter rejection:** `nextflow run . --stage downstream --min_genes 200` fails with an
   unrecognised-parameter error.
4. **Neighbours decoupling (#3):** `qc_clustering` with `--integration_methods harmony` (no scvi)
   either produces clustering output or stops with a clear error.
5. **nf-test** green across the whole suite in CI, on a PR into `dev_ukdri`.
6. **Cluster run:** the three stages on one real UK DRI dataset complete, and every difference
   from the pre-cleanup results (cell counts per QC stage, cluster assignments, top differential
   genes per contrast) is explained by a specific finding above.

**Stage B acceptance:**

1. `nf-core pipelines lint` passes under the non-nf-core ruleset.
2. `nf-core modules lint` passes, with the four patched nf-core modules showing their registered
   patches and no unregistered differences.
3. Steps 1, 2, 3 and 6 of Stage A pass in the new repo.
4. **Parity with `v0.9.0`:** the new repo and `v0.9.0` on the same real dataset give identical
   cell counts per QC stage, cluster assignments and differential-genes tables. Any difference is a
   porting bug.

---

## Appendix A: `test_scdownstream_images.sh`

Container smoke test for the two `nhecker/*` images. Runs as-is on any x86_64 Linux host with
Apptainer. Use it for A8 and consider committing it (e.g. `tests/containers/`) together with the
A7.5 container-resolution check.

````bash
#!/usr/bin/env bash
# Smoke test for the two scdownstream container images.
# Requires: x86_64 Linux, apptainer (or singularity), internet access, ~10 GB free disk.
# Usage:    bash test_scdownstream_images.sh [workdir]
set -uo pipefail

REPORT_IMG=docker://docker.io/nhecker/scanpy-report:1.11.4-coreinf0.4
PYDESEQ2_IMG=docker://docker.io/nhecker/pydeseq2:0.1

APT=$(command -v apptainer || command -v singularity) || { echo "apptainer/singularity not found"; exit 1; }
mkdir -p "${1:-$PWD/scd_image_test}"
W=$(cd "${1:-$PWD/scd_image_test}" && pwd)
mkdir -p "$W/tmp" "$W/render"
cd "$W"
# /tmp is often a small tmpfs; the image build needs several GB
export APPTAINER_TMPDIR="$W/tmp" SINGULARITY_TMPDIR="$W/tmp"

FAILED=0
pass() { echo "PASS: $1"; }
fail() { echo "FAIL: $1"; FAILED=1; }

echo "== 1/4 Pull images into $W"
for pair in "report.sif $REPORT_IMG" "pydeseq2.sif $PYDESEQ2_IMG"; do
    set -- $pair
    if [ -f "$1" ] || "$APT" pull "$1" "$2"; then pass "pull $2"; else fail "pull $2"; fi
done

echo "== 2/4 Report image: packages"
if "$APT" exec --cleanenv --home "$W/render" report.sif bash -c '
    echo "quarto $(quarto --version)"
    python - <<EOF
import importlib.metadata as md
import scanpy, anndata, liana, gprofiler, papermill, nbclient, ipykernel, IPython, matplotlib, yaml, threadpoolctl
for p in ["scanpy", "anndata", "liana", "gprofiler-official", "papermill", "nbclient", "ipykernel", "matplotlib"]:
    print(f"  {p} {md.version(p)}")
EOF'; then pass "report image packages"; else fail "report image packages"; fi

echo "== 3/4 Report image: quarto render with -P parameters"
cat > "$W/render/smoke.qmd" <<'EOF'
---
title: "scdownstream image smoke test"
jupyter: python3
format:
  html:
    embed-resources: true
---

```{python}
#| tags: [parameters]
FILE = "scdownstream.h5ad"
NTOP = 1
```

```{python}
import numpy as np, anndata as ad, scanpy as sc, matplotlib.pyplot as plt
adata = ad.AnnData(np.random.default_rng(0).poisson(1, (100, 50)).astype("float32"))
sc.pp.normalize_total(adata); sc.pp.log1p(adata); sc.pp.pca(adata, n_comps=5)
print(f"PARAM_OK:{NTOP}:{FILE}")
plt.scatter(adata.obsm["X_pca"][:, 0], adata.obsm["X_pca"][:, 1]); plt.show()
```
EOF
(
    cd "$W/render" && rm -f smoke.html
    "$APT" exec --cleanenv --home "$W/render" \
        --env NUMBA_CACHE_DIR=./tmp/numba,MPLCONFIGDIR=./tmp/matplotlib,XDG_CACHE_HOME=./tmp/matplotlib/cache \
        ../report.sif quarto render smoke.qmd -P NTOP:7 -P FILE:dummy.h5ad --to html
)
if [ -s render/smoke.html ] && grep -q "PARAM_OK:7:dummy.h5ad" render/smoke.html && grep -q "data:image/png" render/smoke.html; then
    pass "quarto render (-P reached the notebook, plot embedded)"
else
    fail "quarto render (see output above; render/smoke.html missing or incomplete)"
fi

echo "== 4/4 pydeseq2 image: synthetic differential expression"
cat > "$W/deseq_smoke.py" <<'EOF'
import numpy as np, pandas as pd, anndata, yaml, pydeseq2
from pydeseq2.dds import DeseqDataSet
from pydeseq2.default_inference import DefaultInference
from pydeseq2.ds import DeseqStats
rng = np.random.default_rng(0)
counts = pd.DataFrame(rng.negative_binomial(5, 0.3, (12, 50)),
                      index=[f"s{i}" for i in range(12)], columns=[f"g{j}" for j in range(50)])
counts.iloc[6:, :5] *= 10
meta = pd.DataFrame({"condition": ["A"] * 6 + ["B"] * 6}, index=counts.index)
inference = DefaultInference(n_cpus=1)
dds = DeseqDataSet(counts=counts, metadata=meta, design="~condition",
                   refit_cooks=True, inference=inference, quiet=True)
dds.deseq2()
stats = DeseqStats(dds, contrast=["condition", "B", "A"], inference=inference, quiet=True)
stats.summary()
top = stats.results_df.sort_values("padj").head(10)
print(top.head(6))
print(f"pydeseq2 {pydeseq2.__version__}, anndata {anndata.__version__}")
spiked = {f"g{j}" for j in range(5)}
assert spiked <= set(top.index), f"spiked genes missing from top 10: {spiked - set(top.index)}"
EOF
if "$APT" exec --cleanenv --home "$W" pydeseq2.sif python "$W/deseq_smoke.py" 2> >(grep -v -i warning >&2); then
    pass "pydeseq2 synthetic DE"
else
    fail "pydeseq2 synthetic DE"
fi

echo
if [ "$FAILED" -eq 0 ]; then echo "ALL CHECKS PASSED"; else echo "SOME CHECKS FAILED"; fi
exit "$FAILED"
````

### Running it on a macOS arm64 machine

*(verified 2026-09-25)* The images are `linux/amd64` only. Quarto (Deno/V8) and numba need a
**full x86_64 VM**. Both lighter options failed:
- qemu *user-mode* emulation inside an arm64 VM: Quarto segfaults, and scrublet/numba crashes.
- A Rosetta VM (`vmType: vz`, `rosetta.enabled`): the Jupyter kernel never answers.

Plain Python (pydeseq2, h5ad writes) works under both, but the reports don't.

```bash
# one-off: Lima needs the x86_64 guest agent on Apple Silicon
brew install lima-additional-guestagents
# full x86_64 system emulation with Apptainer preinstalled (new instance; keep the arm64 one)
limactl start --name apptainer-x86 --vm-type=qemu --arch=x86_64 --tty=false template://apptainer
limactl shell apptainer-x86 -- uname -m            # -> x86_64
# copy the smoke test in and run it; keep the work dir on disk, not the 2 GB /tmp tmpfs
limactl copy test_scdownstream_images.sh apptainer-x86:/var/tmp/
limactl shell apptainer-x86 -- bash /var/tmp/test_scdownstream_images.sh /var/tmp/scd_image_test
# clean up
limactl stop apptainer-x86     # or: limactl delete apptainer-x86
```

Result: `ALL CHECKS PASSED` (pulls, report image packages, `quarto render -P` with parameters
reaching the notebook and the plot embedded, and a pydeseq2 DE run). It took **~60 min** wall time
(4 vCPU, 4 GiB), mostly pulling and converting the two 3.3 GB images under emulation; later runs
reuse the `.sif` files. Use `--cpus`/`--memory` on `limactl start` for a faster VM. For running
the pipeline itself, a real x86_64 Linux host is still the practical option.

## Appendix B: container hosting and upstreaming

Public Docker Hub is allowed for this repo's modules (decision 2026-09-01, see `PLAN_DONE.md`).
Seqera Containers is preferred wherever a plain `environment.yml` suffices, and **required** for
any module that is to go to `nf-core/modules`: nf-core's hierarchy is Bioconda + Seqera Containers
first, custom images only under the nf-core quay.io account, other registries only where
licensing restricts redistribution. So:

- `decoupler/pseudobulk` and `decoupler/filter_pseudobulk` use `gcfntnu/scanpy:1.11.4`, a
  third-party image; they need a Seqera Containers image before contributing.
- `pydeseq2/differential_genes` uses the personal `docker.io/nhecker/pydeseq2:0.1`; same.
- The report modules stay local, so Docker Hub is fine for them.

Docker Hub applies anonymous pull rate limits. A shared `NXF_SINGULARITY_CACHEDIR` on the cluster
makes this a non-issue there, but CI runners pulling cold may need an authenticated pull.
