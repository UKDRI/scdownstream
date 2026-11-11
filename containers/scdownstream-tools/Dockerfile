# Image: docker.io/nhecker/scdownstream-tools:1.11.4-1
#
# Build and push (from the repository root):
#   docker build --platform linux/amd64 -t docker.io/nhecker/scdownstream-tools:1.11.4-1 containers/scdownstream-tools && docker push docker.io/nhecker/scdownstream-tools:1.11.4-1
#
# The pipeline's scanpy/anndata generation (docker.io/gcfntnu/scanpy:1.11.4) plus Quarto and the Python
# tools of the h5ad steps that the base does not ship. See containers/README.md for the modules using it.
FROM docker.io/gcfntnu/scanpy:1.11.4

RUN apt update

RUN wget -qO /tmp/quarto-1.9.38-linux-amd64.deb https://github.com/quarto-dev/quarto-cli/releases/download/v1.9.38/quarto-1.9.38-linux-amd64.deb && \
    apt install -y /tmp/quarto-1.9.38-linux-amd64.deb && \
    rm /tmp/quarto-1.9.38-linux-amd64.deb && \
    rm -rf /var/lib/apt/lists/*

# Every package of the base is a constraint, so the tools and their dependencies resolve against the base
# versions and none of them changes. `pip list --format=freeze`, not `pip freeze`: pip freeze lists the conda
# packages as `@ file://` URLs, which pip rejects as constraints.
RUN pip list --format=freeze > /opt/base-constraints.txt && \
    pip install --no-cache-dir -c /opt/base-constraints.txt \
        celltypist==1.6.3 \
        harmonypy==0.0.10 \
        upsetplot==0.9.0 \
        mygene==3.2.2 \
        hugo-unifier==0.3.2 \
        cosg==1.2.1 \
        liana==1.6.1 \
        pydeseq2==0.5.4 \
        papermill==2.7.0 \
        nbconvert==7.17.1 \
        gprofiler-official==1.0.0

RUN python -c "import celltypist, harmonypy, upsetplot, mygene, cosg, liana, pydeseq2, papermill, gprofiler, anndata; assert anndata.__version__ == '0.12.2', anndata.__version__" && \
    hugo-unifier --version && \
    quarto --version
