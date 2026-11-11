# Image: docker.io/nhecker/scdownstream-scimilarity:0.3.0-1
#
# Build and push (from the repository root):
#   docker build --platform linux/amd64 -t docker.io/nhecker/scdownstream-scimilarity:0.3.0-1 containers/scdownstream-scimilarity && docker push docker.io/nhecker/scdownstream-scimilarity:0.3.0-1
#
# SCimilarity 0.3.0 on the pipeline's scanpy/anndata generation (docker.io/gcfntnu/scanpy:1.11.4). torch is the
# default PyPI build with the CUDA 12 libraries, so use_gpu=True (gpu profile) works; it runs on CPU as well.
# See containers/README.md for the modules using it.
FROM docker.io/gcfntnu/scanpy:1.11.4

# hnswlib 0.8.0 has no wheels on PyPI (sdist only), conda-forge has builds. Installed before the pip step, so it
# is part of the constraints there.
RUN mamba install -y --freeze-installed -c conda-forge hnswlib=0.8.0 && \
    mamba clean -a -y

# Every other package of the base is a constraint, so scimilarity's dependencies resolve against the base
# versions. `pip list --format=freeze`, not `pip freeze`: pip freeze lists the conda packages as `@ file://`
# URLs, which pip rejects as constraints.
# zarr and numcodecs are left out of the constraints, because they have to change: scimilarity 0.3.0's
# CellEmbedding.get_embeddings() tests `isinstance(X.data, zarr.core.Array)` for every sparse X, and zarr 3 has
# no zarr.core.Array (AttributeError). anndata 0.12.2 allows zarr >= 2.18.7, != 3.0.*, so zarr 2.18.7 is used;
# it needs numcodecs < 0.16.
RUN pip list --format=freeze | grep -viE '^(zarr|numcodecs)==' > /opt/base-constraints.txt && \
    pip install --no-cache-dir -c /opt/base-constraints.txt \
        torch==2.8.0 \
        scimilarity==0.3.0 \
        captum==0.9.0 \
        pytorch-lightning==2.6.6 \
        tiledb==0.36.1 \
        circlify==0.16.0 \
        obonet==1.3.0 \
        zarr==2.18.7 \
        numcodecs==0.15.1

RUN python -c "import anndata, hnswlib, zarr, scimilarity; from scimilarity import CellAnnotation, CellQuery; from scimilarity.utils import align_dataset, lognorm_counts, pseudobulk_anndata; assert anndata.__version__ == '0.12.2', anndata.__version__; assert zarr.__version__.startswith('2.'), zarr.__version__; print('scimilarity', scimilarity.__version__)"
