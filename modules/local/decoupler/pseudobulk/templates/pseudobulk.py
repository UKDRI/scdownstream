#!/usr/bin/env python3

import os
import platform
import warnings
import yaml

warnings.filterwarnings("ignore")

import anndata as ad
import pandas as pd
import decoupler as dc

adata = ad.read_h5ad("${h5ad}")
prefix = "${prefix}"
sample_col = "${sample_col}"
group_cols = "${group_cols}"

pb_col = "pseudobulk_label"
separator = "___"  

list_group_cols = group_cols.split(',')
required_cols = [sample_col] + list_group_cols

missing = [col for col in required_cols if col not in adata.obs.columns]
if missing:
    raise ValueError(f"Missing required obs columns: {missing}")

adata.obs[pb_col] = adata.obs[list_group_cols].astype(str).agg(separator.join, axis=1)


# decoupler's pseudobulk helper expects a single grouping key.
# We pass the full list and let the function build the grouped object.
pseudobulk_adata = dc.pp.pseudobulk(
    adata,
    sample_col=sample_col,
    groups_col=pb_col,
    layer='counts',
    mode="sum"
)

# re-assign label with separator
pseudobulk_adata.obs_names = pseudobulk_adata.obs[sample_col].astype(str) + separator + pseudobulk_adata.obs[pb_col].astype(str)


pseudobulk_adata.write_h5ad(f"{prefix}.h5ad")

versions = {
    "${task.process}": {
        "python": platform.python_version(),
        "anndata": ad.__version__,
        "decoupler": dc.__version__,
        "pandas": pd.__version__,
    }
}
with open("versions.yml", "w") as f:
    yaml.dump(versions, f)
