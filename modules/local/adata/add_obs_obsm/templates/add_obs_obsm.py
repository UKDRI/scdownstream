#!/usr/bin/env python3

import platform
import os

import anndata as ad
import pandas as pd
import numpy as np
import yaml

adata = ad.read_h5ad("${adata_base}")
adata_add = ad.read_h5ad("${adata_add}")
prefix = "${prefix}"


for col in adata_add.obs:
    if col not in adata.obs:
        adata.obs[col] = adata_add.obs[col].copy()

for emb in adata_add.obsm:
    if emb not in adata.obsm:
        adata.obsm[emb] = adata_add.obsm[emb].copy()


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