# Image: docker.io/nhecker/scanpy-cosg:1.11.4-cosg1.2.1
#
FROM docker.io/gcfntnu/scanpy:1.11.4

# --no-deps: COSG's dependencies (anndata, numpy, pandas, scipy, scikit-learn, numba, networkx, matplotlib) come from the base image, so the versions stay those of the rest of the pipeline
RUN pip install --no-cache-dir --no-deps cosg==1.2.1 && python -c "import cosg"
