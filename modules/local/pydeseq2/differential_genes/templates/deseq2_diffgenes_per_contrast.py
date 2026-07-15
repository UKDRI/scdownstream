#!/usr/bin/env python3

import platform
import warnings
import yaml

warnings.filterwarnings("ignore")

import anndata as ad
import pandas as pd
import decoupler as dc
from pydeseq2.dds import DeseqDataSet
from pydeseq2.ds import DeseqStats

adata = ad.read_h5ad("${h5ad}")
prefix = "${prefix}"
condition_col = "${condition_col}"
contrast = "${contrast}"

if condition_col not in adata.obs.columns:
    raise ValueError(f"Condition column '{condition_col}' not found in adata.obs")

# Use pseudobulk AnnData directly and compare the requested contrast.
contrast_mask = adata.obs[condition_col].astype(str) == contrast
if contrast_mask.sum() < 3:
    raise ValueError(f"Contrast '{contrast}' has too few samples ({contrast_mask.sum()})")

other_mask = ~contrast_mask
if other_mask.sum() < 3:
    raise ValueError(f"Reference group for contrast '{contrast}' has too few samples ({other_mask.sum()})")

counts = pd.DataFrame(
    adata.X.toarray() if hasattr(adata.X, "toarray") else adata.X,
    index=adata.var_names,
    columns=adata.obs_names,
)
if hasattr(adata.layers, "get") and "counts" in adata.layers:
    counts = pd.DataFrame(
        adata.layers["counts"].toarray() if hasattr(adata.layers["counts"], "toarray") else adata.layers["counts"],
        index=adata.var_names,
        columns=adata.obs_names,
    )

counts = counts.loc[:, contrast_mask | other_mask]
counts = counts.apply(pd.to_numeric, errors="coerce").fillna(0)
counts = counts.clip(lower=0).round().astype(int)

metadata = pd.DataFrame(index=counts.columns)
metadata[condition_col] = [contrast if flag else "reference" for flag in (contrast_mask.loc[counts.columns].values)]
metadata["contrast_group"] = ["contrast" if flag else "reference" for flag in (contrast_mask.loc[counts.columns].values)]

dds = DeseqDataSet(
    counts=counts.T,
    metadata=metadata,
    design_factors=["contrast_group"],
    refit_cooks=True,
    n_cpus=1,
)
dds.deseq2()
dds.normalize()
stats = DeseqStats(dds, alpha=0.05)
stats.summary()

result_df = stats.results_df.reset_index().rename(columns={"index": "feature"})
result_df.to_csv(f"{prefix}.csv", index=False)

versions = {
    "${task.process}": {
        "python": platform.python_version(),
        "anndata": ad.__version__,
        "decoupler": dc.__version__,
        "pydeseq2": __import__("pydeseq2").__version__,
        "pandas": pd.__version__,
    }
}
with open("versions.yml", "w") as f:
    yaml.dump(versions, f)
