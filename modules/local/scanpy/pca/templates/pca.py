#!/usr/bin/env python3

import os
import platform

os.environ["NUMBA_CACHE_DIR"] = "./tmp/numba"
os.environ["MPLCONFIGDIR"] = "./tmp/matplotlib"

import scanpy as sc
import numpy as np
import pandas as pd
import yaml

from threadpoolctl import threadpool_limits
threadpool_limits(int("${task.cpus}"))
sc.settings.n_jobs = int("${task.cpus}")

adata = sc.read_h5ad("${h5ad}")
prefix = "${prefix}"
n_comps = "${task.ext.n_comps ?: ''}"
# Fixed seed keeps the PCA reproducible between runs
random_state = 0
# Rounding the PCA coordinates keeps output hashes stable
n_decimals = 8

if n_comps:
    n_comps = int(n_comps)
    # scanpy runs PCA on the highly variable genes when var["highly_variable"] is set
    n_genes = int(adata.var["highly_variable"].sum()) if "highly_variable" in adata.var else adata.n_vars
    n_max = min(adata.n_obs, n_genes) - 1
    if n_comps > n_max:
        print(f"WARNING: pca_n_comps={n_comps} exceeds what this object allows (min(n_obs, genes used for PCA) - 1 = {n_max}; genes used for PCA are the highly variable genes when set), using {n_max}")
        n_comps = n_max


# Run PCA
sc.pp.pca(adata, n_comps=n_comps or None, random_state=random_state)

# Round to n_decimals decimal places
# This ensures hashes are stable
key_added = 'X_pca'
adata.obsm[key_added] = np.round(adata.obsm[key_added], n_decimals)

adata.write_h5ad(f"{prefix}.h5ad")
df = pd.DataFrame(adata.obsm[key_added], index=adata.obs_names)
df.to_pickle(f"X_{prefix}.pkl")

# Versions
versions = {
    "${task.process}": {
        "python": platform.python_version(),
        "scanpy": sc.__version__,
        "pandas": pd.__version__
    }
}

with open("versions.yml", "w") as f:
    yaml.dump(versions, f)
