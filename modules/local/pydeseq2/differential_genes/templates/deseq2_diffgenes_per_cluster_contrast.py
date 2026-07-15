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

clusters = sorted([str(v) for v in adata.obs[cluster_col].dropna().astype(str).unique().tolist()])

for cluster_value in clusters:
    cluster_mask = adata.obs[cluster_col].astype(str) == cluster_value
    if cluster_mask.sum() < 3:
        print(f"Skipping cluster '{cluster_value}' because it has too few samples ({cluster_mask.sum()})")
        continue

    for contrast_value in sorted([str(v) for v in adata.obs[condition_col].dropna().astype(str).unique().tolist()]):
        contrast_mask = adata.obs[condition_col].astype(str) == contrast_value
        if contrast_mask.sum() < 3:
            continue

        group_mask = cluster_mask & contrast_mask
        if group_mask.sum() < 3:
            continue

        other_mask = (~cluster_mask) & (~contrast_mask)
        if other_mask.sum() < 3:
            continue

        sample_mask = group_mask | other_mask
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

        counts = counts.loc[:, sample_mask]
        counts = counts.apply(pd.to_numeric, errors="coerce").fillna(0)
        counts = counts.clip(lower=0).round().astype(int)

        metadata = pd.DataFrame(index=counts.columns)
        metadata[condition_col] = [contrast_value if flag else "reference" for flag in (sample_mask.loc[counts.columns].values)]
        metadata["contrast_group"] = ["cluster_contrast" if flag else "reference" for flag in (sample_mask.loc[counts.columns].values)]

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
        safe_contrast = str(contrast_value).replace(os.sep, "_").replace(" ", "_")
        result_df.to_csv(f"{prefix}_{safe_cluster}_{safe_contrast}.csv", index=False)

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
