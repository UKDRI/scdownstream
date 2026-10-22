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
exclude_col = "${exclude_samples_col}"
exclude_values_raw = "${exclude_samples_values}"


def sanitize(name):
    """Replace characters that break the formula / file naming with underscores."""
    name = str(name)
    for c in ['~', ':', '+', '-', ' ']:
        name = name.replace(c, '_')
    return name


# Per-contrast sample exclusion (nf-core/differentialabundance semantics): for this
# contrast only, drop every pseudobulk sample whose value in `exclude_col` is one of the
# comma-separated `exclude_values`. Applied before the counts, metadata and replicate
# guard are built, so every later step (including the min_samples check) sees only the
# kept samples.
exclude_col = exclude_col.strip()
if exclude_col == "null":
    exclude_col = ""
exclude_values = []
if exclude_values_raw and exclude_values_raw.strip() not in ("", "null"):
    exclude_values = [v.strip() for v in exclude_values_raw.split(',') if v.strip() != ""]

n_excluded = 0
if exclude_col and not exclude_values:
    print(
        f"WARNING: contrast '{contrast_name}': exclude_samples_col '{exclude_col}' is set "
        f"but exclude_samples_values is empty; no samples excluded"
    )
elif exclude_values and not exclude_col:
    print(
        f"WARNING: contrast '{contrast_name}': exclude_samples_values {exclude_values} are set "
        f"but exclude_samples_col is empty; no samples excluded"
    )
elif exclude_col and exclude_values:
    if exclude_col not in adata.obs.columns:
        # The object is the pseudobulk profile set (one row per sample x group label).
        # Pseudobulk aggregation keeps only obs columns that are constant within each
        # profile, so columns that vary within a sample are gone at this point.
        raise ValueError(
            f"Contrast '{contrast_name}': exclude_samples_col '{exclude_col}' not found in "
            f"the pseudobulk obs columns. Only obs columns that are constant within each "
            f"pseudobulk profile survive aggregation, so this must be sample-level metadata "
            f"(e.g. donor, batch, sex). Available columns: "
            f"{', '.join(map(str, adata.obs.columns))}"
        )
    exclude_series = adata.obs[exclude_col].astype(str)
    mask = exclude_series.isin(exclude_values).to_numpy()
    present_values = set(exclude_series)
    unmatched = [v for v in exclude_values if v not in present_values]
    n_excluded = int(mask.sum())
    adata = adata[~mask].copy()
    print(
        f"Contrast '{contrast_name}': excluding samples with {exclude_col} in "
        f"[{', '.join(exclude_values)}]: {n_excluded} pseudobulk sample(s) excluded, "
        f"{adata.n_obs} remaining"
    )
    if unmatched:
        print(
            f"WARNING: contrast '{contrast_name}': exclude_samples_values "
            f"[{', '.join(unmatched)}] matched no sample in column '{exclude_col}' "
            f"(possible typo)"
        )


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
    # Name the exclusion, so a skip caused by it is understandable from the log.
    exclusion_note = ""
    if n_excluded > 0:
        exclusion_note = (
            f" after excluding {n_excluded} sample(s) with {exclude_col} in "
            f"[{', '.join(exclude_values)}]"
        )
    print(
        f"Skipping contrast '{contrast_name}': "
        f"{describe('target', target_group, n_target)}, "
        f"{describe('reference', reference_group, n_reference)} "
        f"(minimum {min_samples} per group){exclusion_note}"
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
