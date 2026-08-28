#!/usr/bin/env python3
import os
import platform
import warnings
import yaml

warnings.filterwarnings("ignore")

import anndata as ad
import pandas as pd
from pydeseq2.dds import DeseqDataSet
from pydeseq2.ds import DeseqStats
from pydeseq2.default_inference import DefaultInference

adata = ad.read_h5ad("${h5ad}")
prefix = "${prefix}"
ncpus = int("${task.cpus}")

cluster_col = "${cluster_col}"
blocked_raw = "${blocking_col}"


def sanitize(name):
    """Replace characters that break the formula / file naming with underscores."""
    name = str(name)
    for c in ['~', ':', '+', '-', ' ']:
        name = name.replace(c, '_')
    return name


if cluster_col not in adata.obs.columns:
    raise ValueError(f"Cluster column '{cluster_col}' not found in adata.obs")

# Build the counts matrix once (samples x genes, integer counts expected by pydeseq2).
counts = pd.DataFrame(
    adata.X.toarray() if hasattr(adata.X, "toarray") else adata.X,
    index=adata.obs_names,
    columns=adata.var_names,
)
counts = counts.clip(lower=0).round().astype(int)

# Optional blocking covariate(s), added to the design formula and checked for existence.
list_blocked = []
block_metadata = pd.DataFrame(index=counts.index)
if blocked_raw and blocked_raw.strip() not in ("", "null"):
    for col in blocked_raw.split(','):
        col = col.strip()
        if col == "":
            continue
        if col not in adata.obs.columns:
            raise ValueError(f"Blocking column '{col}' not found in adata.obs")
        colname = sanitize(col)
        block_metadata[colname] = adata.obs[col].values
        list_blocked.append(colname)

# Sanitized, formula-safe name for the one-vs-rest membership factor.
factor = "cluster_membership"

# One-vs-rest: test each cluster level against all other pooled samples.
cluster_series = adata.obs[cluster_col].astype(str)
levels = sorted(cluster_series.dropna().unique().tolist())

for level in levels:
    this_mask = cluster_series == level
    n_this = int(this_mask.sum())
    n_rest = int((~this_mask).sum())

    if n_this < 3:
        print(f"Skipping cluster '{level}': too few samples in group ({n_this} < 3)")
        continue
    if n_rest < 3:
        print(f"Skipping cluster '{level}': too few samples in rest ({n_rest} < 3)")
        continue

    # Two-level design factor: "this" for group members, "rest" for everyone else.
    metadata = block_metadata.copy()
    metadata[factor] = ["this" if flag else "rest" for flag in this_mask]

    # Compile design formula: ~ cluster_membership [+ blocked1 + blocked2 ...]
    design = "~" + factor
    for blocked in list_blocked:
        design += " + " + blocked

    inference = DefaultInference(n_cpus=ncpus)

    dds = DeseqDataSet(
        counts=counts,
        metadata=metadata,
        design=design,
        refit_cooks=True,
        inference=inference,
    )
    dds.deseq2()

    # Log fold change of this cluster vs the pooled rest.
    stats = DeseqStats(
        dds,
        contrast=[factor, "this", "rest"],
        inference=inference,
    )
    stats.summary()

    result_df = stats.results_df.reset_index().rename(columns={"index": "feature"})
    safe_level = sanitize(level).replace(os.sep, "_")
    result_df.to_csv(f"{prefix}_{safe_level}.tsv", sep="\\t", index=False)

versions = {
    "${task.process}": {
        "python": platform.python_version(),
        "anndata": ad.__version__,
        "pydeseq2": __import__("pydeseq2").__version__,
        "pandas": pd.__version__,
    }
}
with open("versions.yml", "w") as f:
    yaml.dump(versions, f)
