# UK DRI scdownstream

## Use the repo agents for code changes

Implement pipeline changes through the two agents in `.claude/agents/`, not directly. This holds
for small changes too, and for follow-up fixes to code an agent wrote.

- **`pipeline-developer`**: `modules/local/**` (process bodies, Python/R templates,
  `environment.yml`), the Quarto report templates, and the `withName:` entries in
  `conf/modules.config`.
- **`nextflow-dsl2`**: `main.nf`, `workflows/*.nf`, `subworkflows/local/**`, the params in
  `nextflow.config`, `nextflow_schema.json` and `conf/test.config`.

When a change touches both, fix the interface first (process name, input order, emit names,
`task.ext.*` keys) and run the two agents in parallel. Afterwards, review and test their changes,
and update the docs (`README.md`, `docs/usage.md`, `docs/output.md`, `CHANGELOG.md`).

An agent that is resumed with a follow-up message keeps the instructions it started with. If
its definition file changed in the meantime, include the new rule in the message.

## Development gotchas

- **One anndata/pandas generation for every container.** Processes hand h5ad files to each other,
  so every image that reads or writes an h5ad has to use the pipeline's versions: the base image
  `docker.io/gcfntnu/scanpy:1.11.4` (python 3.11, scanpy 1.11.4, anndata 0.12.2, pandas 2.3.x,
  numpy 2.2). **pandas ≥ 3 / anndata ≥ 0.13 are not compatible.** Older anndata (0.10, 0.11) reads
  most of what 0.12 writes, but not `None` values: anndata ≥ 0.12 stores them as `null` elements.
  - *What happened (2026-11-10):* `MARKERS_COSG` ran in the BioContainer `cosg:1.2.1`. It is built
    with unpinned dependencies and so had pandas 3.1 / anndata 0.13. Those write string indices as
    nullable string arrays. `SCANPY_ENRICH` (anndata 0.12.2) then failed:
    `allow_write_nullable_strings is False ... /obs`.
  - *What happened (2026-11-11):* the presto/COSG templates stored `"layer": None` in
    `uns["rank_genes_groups"]["params"]`. `LIANA_RANKAGGREGATE` (BioContainer `liana:1.4.0`, python
    3.9, anndata 0.10) then failed: `No read method registered for IOSpec(encoding_type='null'...)`.
    `ADATA_TORDS` (anndata 0.11.1) and `ADATA_SPLITCOL` (0.11.4) failed the same way.
  - *Rule:* run every h5ad-handling process on the base image or an image built `FROM` it. scanpy
    itself writes `None` (e.g. `sc.pp.log1p` stores `uns["log1p"]["base"] = None`), so avoiding
    `None` in our templates is not enough. Since 2026-11-11 (`anndata-0.12-throughout`) the main
    path is on the base image and `ADATA_TORDS` is removed; the modules still on older images are
    listed in PLAN.md. Don't store `None` in an h5ad in our templates anyway; leave the key out.
  - *Rule:* don't use BioContainers / Seqera images with unpinned or newer anndata/pandas for
    h5ad-handling processes. Extend `gcfntnu/scanpy:1.11.4` instead, with
    `pip install --no-deps <tool>` (see `modules/local/markers/cosg/Dockerfile`), or pin the versions.
  - *Upgrading:* moving to pandas 3 / anndata 0.13 means updating all images together, then testing
    every stage end to end.
