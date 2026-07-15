#!/usr/bin/env python3

import os
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
cluster_col = "${cluster_col}"
condition_col = "${condition_col}"

if cluster_col not in adata.obs.columns:
    raise ValueError(f"Cluster column '{cluster_col}' not found in adata.obs")
if condition_col not in adata.obs.columns:
    raise ValueError(f"Condition column '{condition_col}' not found in adata.obs")

# Use decoupler to prepare pseudobulk-compatible metadata and keep the input as pseudobulk AnnData.
clusters = sorted([str(v) for v in adata.obs[cluster_col].dropna().astype(str).unique().tolist()])

for cluster_value in clusters:
    cluster_mask = adata.obs[cluster_col].astype(str) == cluster_value
    if cluster_mask.sum() < 3:
        print(f"Skipping cluster '{cluster_value}' because it has too few samples ({cluster_mask.sum()})")
        continue

    sub_adata = adata[cluster_mask].copy()
    if sub_adata.n_obs < 3:
        continue

    # Compare this cluster against all other pseudobulk samples.
    other_mask = ~cluster_mask
    other_sub_adata = adata[other_mask].copy()
    if other_sub_adata.n_obs < 3:
        continue

    contrast_values = pd.Series(["cluster" if flag else "other" for flag in cluster_mask], index=adata.obs_names)
    contrast_values = contrast_values[cluster_mask | other_mask]

    counts = pd.DataFrame(
        sub_adata.X.toarray() if hasattr(sub_adata.X, "toarray") else sub_adata.X,
        index=sub_adata.var_names,
        columns=sub_adata.obs_names,
    )
    other_counts = pd.DataFrame(
        other_sub_adata.X.toarray() if hasattr(other_sub_adata.X, "toarray") else other_sub_adata.X,
        index=other_sub_adata.var_names,
        columns=other_sub_adata.obs_names,
    )

    counts = pd.concat([counts, other_counts], axis=1)
    counts = counts.apply(pd.to_numeric, errors="coerce").fillna(0)
    counts = counts.clip(lower=0).round().astype(int)

    metadata = pd.DataFrame(index=counts.columns)
    metadata[condition_col] = [cluster_value, "other"] if len(counts.columns) == 2 else [cluster_value] + ["other"] * (counts.shape[1] - 1)
    metadata["contrast_group"] = ["cluster", "other"] if len(counts.columns) == 2 else ["cluster"] + ["other"] * (counts.shape[1] - 1)

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
    safe_cluster = str(cluster_value).replace(os.sep, "_").replace(" ", "_")
    result_df.to_csv(f"{prefix}_{safe_cluster}.csv", index=False)

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
