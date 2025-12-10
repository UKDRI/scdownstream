#!/usr/bin/env python3

import platform
import os

import anndata as ad
import pandas as pd
import numpy as np
import yaml

# TODO:  extract obs, obsm, obsnanes for adata_add, del adata_add, get index for subsetting from shared genes
adata = ad.read_h5ad("${adata_base}")
adata_add = ad.read_h5ad("${adata_add}")
prefix = "${prefix}"

# get obs, obsm
add_obs = adata_add.obs.copy()
add_obsm = adata_add.obsm.copy()
add_obsnames = [str(bc) for bc in adata_add.obs_names ]
del adata_add

# subset to overlapping cells
shared_bcs = [bc for bc in adata.obs_names if bc in add_obsnames ]

# subset to shared and add to base adata
adata = adata[shared_bcs].copy()

dict_idx = {bc:idx for idx, bc in enumerate(add_obsnames) }
order_add = [ dict_idx[bc] for bc in shared_bcs ]
add_obs = add_obs.iloc[order_add, :].copy()

for col in add_obs:
    if col not in adata.obs:
        adata.obs[col] = add_obs[col].copy()
for emb in add_obsm:
    if emb not in adata.obsm:
        adata.obsm[emb] = add_obsm[emb][order_add, :].copy()


adata.write_h5ad(f"{prefix}.h5ad")

# Versions

versions = {
    "${task.process}": {
        "python": platform.python_version(),
        "anndata": ad.__version__,
        "pandas": pd.__version__,
        "numpy": np.__version__
    }
}

with open("versions.yml", "w") as f:
    yaml.dump(versions, f)