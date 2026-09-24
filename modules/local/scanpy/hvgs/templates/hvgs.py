#!/usr/bin/env python3

import os
import platform
from threadpoolctl import threadpool_limits

os.environ["MPLCONFIGDIR"] = "./tmp/mpl"
os.environ["NUMBA_CACHE_DIR"] = "./tmp/numba"

import scanpy as sc
import yaml

threadpool_limits(int("${task.cpus}"))
sc.settings.n_jobs = int("${task.cpus}")

adata = sc.read_h5ad("${h5ad}")
prefix = "${prefix}"
n_hvgs = int("${n_hvgs}")
batch_key = "${batch_key}"
subset_to_hvgs = "${subset_to_hvgs}" == "true"


select_hvgs = adata.n_vars > n_hvgs and n_hvgs >= 0

kwargs = {}

if batch_key:
    kwargs["batch_key"] = batch_key

# If an actual limit is provided, use it
# Otherwise, scanpy will automatically determine the number of highly variable genes
if select_hvgs and n_hvgs > 0:
    kwargs["n_top_genes"] = n_hvgs

# Always run this: it also stores the statistics (uns["hvg"], means, dispersions) that
# sc.pl.highly_variable_genes needs in the reports
sc.pp.highly_variable_genes(adata, **kwargs)

if not select_hvgs:
    # Not more genes than requested (or selection disabled): keep every gene
    adata.var["highly_variable"] = True

adata.var[["highly_variable"]].to_pickle(f"{prefix}.pkl")

if subset_to_hvgs:
    adata = adata[:, adata.var["highly_variable"]]

adata.write_h5ad(f"{prefix}.h5ad")

# Versions

versions = {
    "${task.process}": {
        "python": platform.python_version(),
        "scanpy": sc.__version__
    }
}

with open("versions.yml", "w") as f:
    yaml.dump(versions, f)
