# Container images

The processes hand h5ad files to each other, so every image that reads or writes an h5ad uses the pipeline's
anndata/pandas generation: the base image `docker.io/gcfntnu/scanpy:1.11.4` (python 3.11, scanpy 1.11.4,
anndata 0.12.2, pandas 2.3.2, numpy 2.2.6). The images here add the tools the base does not ship.

| Image | Dockerfile | Modules |
| --- | --- | --- |
| `docker.io/nhecker/scanpy-ukdri-suite:1.11.4-1` | `scanpy-ukdri-suite/Dockerfile` | `adata/mygene`, `adata/upsetgenes`, `celltypes/celltypist`, `celltypes/celltypist_train`, `doublet_detection/doublet_removal`, `hugounifier/apply`, `hugounifier/get`, `liana/rankaggregate`, `markers/cosg`, `pydeseq2/differential_genes`, `pydeseq2/report`, `scanpy/enrich`, `scanpy/harmony`, `scanpy/report` |
| `docker.io/nhecker/seurat-ukdri-suite:5.5.1-1` | `seurat-ukdri-suite/Dockerfile` (+ `environment.yml`) | `celda/decontx`, `doublet_detection/scds`, `seurat/integration`, `soupx` |

The image names are deliberately not tied to the pipeline or repository name, so they stay valid if the
pipeline moves. Tags are `<main tool version>-<build number>` (scanpy for `scanpy-ukdri-suite`, Seurat for
`seurat-ukdri-suite`); the build number increases when the image is rebuilt with the same main tool version.

`scanpy-ukdri-suite` contains: Quarto 1.9.38, celltypist 1.6.3, harmonypy 0.0.10, upsetplot 0.9.0,
mygene 3.2.2, hugo-unifier 0.3.2, cosg 1.2.1, liana 1.6.1, pydeseq2 0.5.4, papermill 2.7.0, nbconvert 7.17.1
and gprofiler-official 1.0.0.

`seurat-ukdri-suite` is not built on the base: anndataR 1.2.1 needs R >= 4.5 and rhdf5 >= 2.52.1, and the base has
R 4.4.1. It is a micromamba image like the Wave images, with R 4.5.3, anndataR 1.2.1, celda 1.26.0,
scds 1.26.1, SoupX 1.6.2 and Seurat 5.5.1. anndataR is pinned to exactly 1.2.1: 1.2.2 and later write R
character vectors containing `NA` as `nullable-string-array`, which anndata 0.12.2 reads but refuses to write
back (`allow_write_nullable_strings is False`).

## Build and push

From the repository root, for each image (example for `scanpy-ukdri-suite`):

```bash
docker build --platform linux/amd64 -t docker.io/nhecker/scanpy-ukdri-suite:1.11.4-1 containers/scanpy-ukdri-suite
docker push docker.io/nhecker/scanpy-ukdri-suite:1.11.4-1
```

The tag of each image is in the `# Image:` line at the top of its Dockerfile. A changed image gets a new
tag (a new main tool version, or the same one with the build number increased), and the modules are
pointed at it.

## Rules for a new or changed image

- Build `FROM docker.io/gcfntnu/scanpy:1.11.4` for anything that reads or writes h5ad. Don't use
  BioContainers or Seqera/Wave images with unpinned or newer anndata/pandas (pandas >= 3, anndata >= 0.13).
- Pin the version of every tool that is installed.
- Install Python packages with the base packages as constraints, so that no base package changes:

  ```dockerfile
  RUN pip list --format=freeze > /opt/base-constraints.txt && \
      pip install --no-cache-dir -c /opt/base-constraints.txt <tool>==<version>
  ```

  Use `pip list --format=freeze`, not `pip freeze`: pip freeze lists the conda packages as `@ file://` URLs,
  which pip rejects as constraints. Leave a package out of the constraints only when it has to change, and
  say why in a comment.
- End with a check step that imports the tools and asserts the anndata version.
- Update the `environment.yml` of the modules using the image to the image's versions.
