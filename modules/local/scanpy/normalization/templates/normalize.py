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

# Saving count data in layer or restoring it to avoid double normalization
if 'counts' not in adata.layers:
    adata.layers["counts"] = adata.X.copy()
else:
    adata.X = adata.layers["counts"].copy()

# Normalizing to median total counts
sc.pp.normalize_total(adata)
# Logarithmize the data
sc.pp.log1p(adata)

# write output
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