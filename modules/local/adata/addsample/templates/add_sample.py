#!/usr/bin/env python3

import platform

import anndata as ad
import pandas as pd
import yaml

# The input is a symlink into the previous task's work directory; never write through it
if "${h5ad}" == "${prefix}.h5ad":
    raise ValueError("ADATA_ADDSAMPLE: input and output are both named ${prefix}.h5ad; set a different ext.prefix")

sample_id = "${meta.id}"

adata = ad.read_h5ad("${h5ad}")

# Keep a differing pre-existing "sample" column as "sample_original"
if "sample" in adata.obs:
    old_sample = adata.obs["sample"].astype(str)
    if (old_sample == sample_id).all():
        print(f"obs['sample'] already equals '{sample_id}' for all cells; leaving it as is")
    else:
        if "sample_original" in adata.obs:
            raise ValueError(
                f"ADATA_ADDSAMPLE: obs['sample'] holds values other than '{sample_id}', "
                "but obs['sample_original'] already exists, so it cannot be kept. "
                "Rename or drop one of the two columns in the input h5ad."
            )
        distinct = list(old_sample.unique())
        shown = ", ".join(distinct[:10]) + (f", ... ({len(distinct)} in total)" if len(distinct) > 10 else "")
        print(f"WARNING: obs['sample'] held other values ({shown}); moved to obs['sample_original'] and set obs['sample'] to '{sample_id}'")
        adata.obs["sample_original"] = adata.obs["sample"]

adata.obs["sample"] = pd.Categorical([sample_id] * adata.n_obs)

adata.write_h5ad("${prefix}.h5ad")

# Versions

versions = {
    "${task.process}": {
        "python": platform.python_version(),
        "anndata": ad.__version__,
        "pandas": pd.__version__,
    }
}

with open("versions.yml", "w") as f:
    yaml.dump(versions, f)
