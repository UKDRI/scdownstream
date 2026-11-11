# Image: docker.io/nhecker/scdownstream-r:4.5-1
#
# Build and push (from the repository root):
#   docker build --platform linux/amd64 -t docker.io/nhecker/scdownstream-r:4.5-1 containers/scdownstream-r && docker push docker.io/nhecker/scdownstream-r:4.5-1
#
# R tools of the h5ad steps (anndataR, celda, scds, SoupX, Seurat). anndataR 1.2.1 needs R >= 4.5 and
# rhdf5 >= 2.52.1, so it cannot go into the R 4.4.1 of docker.io/gcfntnu/scanpy:1.11.4. The layout follows the
# Wave images (same micromamba base). See containers/README.md for the modules using it.
FROM mambaorg/micromamba:1.5.10-noble

COPY --chown=$MAMBA_USER:$MAMBA_USER environment.yml /tmp/environment.yml

RUN micromamba install -y -n base -f /tmp/environment.yml && \
    micromamba clean -a -y

USER root
ENV PATH="$MAMBA_ROOT_PREFIX/bin:$PATH"

RUN Rscript -e 'library(anndataR); library(rhdf5); library(celda); library(scds); library(SoupX); library(Seurat); library(glmGamPoi); stopifnot(packageVersion("anndataR") == "1.2.1")'
