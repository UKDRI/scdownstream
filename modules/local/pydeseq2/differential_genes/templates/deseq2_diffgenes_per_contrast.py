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

# Minimum number of pseudobulk samples required in *each* level of the contrast.
min_samples = int("${min_samples}")

# Single contrast definition, supplied one row per task by the caller.
contrast_name = "${contrast_name}"
variable = "${variable}"
target_group = "${target_group}"
reference_group = "${reference_group}"
blocked_raw = "${blocked_variables}"


def sanitize(name):
    """Replace characters that break the formula / file naming with underscores."""
    name = str(name)
    for c in ['~', ':', '+', '-', ' ']:
        name = name.replace(c, '_')
    return name


# Build the counts matrix (samples x genes, integer counts expected by pydeseq2).
counts = pd.DataFrame(
    adata.X.toarray() if hasattr(adata.X, "toarray") else adata.X,
    index=adata.obs_names,
    columns=adata.var_names,
)

if variable not in adata.obs.columns:
    raise ValueError(f"Variable column '{variable}' not found in adata.obs")

# Sanitized design-factor name, used consistently in metadata, formula and contrast.
var_col = sanitize(variable)

metadata = pd.DataFrame(index=counts.index)
metadata[var_col] = adata.obs[variable].astype(str).values

list_blocked = []
if blocked_raw and blocked_raw.strip() not in ("", "null"):
    for col in blocked_raw.split(','):
        col = col.strip()
        if col == "":
            continue
        if col not in adata.obs.columns:
            raise ValueError(f"Blocked column '{col}' not found in adata.obs")
        colname = sanitize(col)
        metadata[colname] = adata.obs[col].values
        list_blocked.append(colname)

# Replicate guard. The threshold applies *per contrast level*: both the target and the
# reference side must have at least `min_samples` pseudobulk samples on their own. A
# group with plenty of samples overall but only one on one side of the contrast is still
# skipped, because DESeq2 cannot estimate dispersion without replicates on both sides.
# A level that is entirely absent from the column counts as zero samples.
level_counts = metadata[var_col].value_counts()
n_target = int(level_counts.get(target_group, 0))
n_reference = int(level_counts.get(reference_group, 0))


def describe(role, level, n):
    """Human-readable count for one side of the contrast."""
    if level not in level_counts.index:
        return f"{role} '{level}' is absent from column '{variable}' (0 sample(s))"
    return f"{role} '{level}' has {n} sample(s)"


run_analysis = n_target >= min_samples and n_reference >= min_samples

if not run_analysis:
    print(
        f"Skipping contrast '{contrast_name}': "
        f"{describe('target', target_group, n_target)}, "
        f"{describe('reference', reference_group, n_reference)} "
        f"(minimum {min_samples} per group)"
    )
else:
    # Compile design formula: ~ variable [+ blocked1 + blocked2 ...]
    design = "~" + var_col
    for blocked in list_blocked:
        design += " + " + blocked

    # Compute fits.
    inference = DefaultInference(n_cpus=ncpus)

    dds = DeseqDataSet(
        counts=counts,
        metadata=metadata,
        design=design,
        refit_cooks=True,
        inference=inference,
    )
    dds.deseq2()

    # Compute stats for this contrast: log fold change of target vs reference.
    stats = DeseqStats(
        dds,
        contrast=[var_col, target_group, reference_group],
        inference=inference,
    )
    stats.summary()

    result_df = stats.results_df.reset_index().rename(columns={"index": "feature"})
    safe_name = sanitize(contrast_name).replace(os.sep, "_")
    result_df.to_csv(f"{prefix}_{safe_name}.tsv", sep="\\t", index=False)

# versions.yml is a mandatory output, so it is written on both the analysed and the
# skipped path (no early sys.exit above).

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
