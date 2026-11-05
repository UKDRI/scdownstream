#!/usr/bin/env python3

import os
import platform
import sys

os.environ["NUMBA_CACHE_DIR"] = "./tmp/numba"
os.environ["MPLCONFIGDIR"] = "./tmp/matplotlib"

import anndata as ad
import numpy as np
import pandas as pd
import scipy
import scipy.sparse as sp
import yaml

from threadpoolctl import threadpool_limits
threadpool_limits(int("${task.cpus}"))

adata = ad.read_h5ad("${h5ad}")
prefix = "${prefix}"
batch_col = "${batch_col}"
symbol_col = "${symbol_col}"
subset_col = "${subset_col}"
subset_values = "${subset_values}"

# Non-integer values that stay at or below this maximum look log-normalised, not like raw counts
max_lognorm_value = 20

# Names the annotation modules read (CELLTYPES_CELLTYPIST, SCVITOOLS_SCANVI with the batch key of the
# scANVI models); they define the structure of the output, not a setting
counts_key = "counts"
batch_key = "batch"
# Batches when the batch column is missing: the samples, else one batch for all cells
fallback_batch_col = "sample"
single_batch = "batch_1"
# Columns of earlier predictions are dropped, so the annotation modules only return the new ones
prediction_prefixes = ("celltypist:", "scanvi:")
# var column that keeps the input var_names when the gene symbols come from another column
input_var_names_key = "var_names_input"


def warn(message):
    print(f"WARNING: ADATA_PREPANNOTATE: {message}", file=sys.stderr)


def short_list(values, n=20):
    values = [str(v) for v in values]
    return ", ".join(values[:n]) + (f", ... ({len(values)} in total)" if len(values) > n else "")


def print_counts(title, values):
    counts = values.value_counts(sort=False).sort_index()
    print(f"{title}:")
    for value, n in counts.items():
        print(f"  {value}: {n}")


if not adata.obs_names.is_unique:
    raise ValueError(
        "ADATA_PREPANNOTATE: the cell names (obs_names) are not unique; the predictions are added back "
        "to the input by cell name, so they must be unique."
    )

# Cell subset

if subset_col or subset_values:
    if not subset_col or not subset_values:
        raise ValueError("ADATA_PREPANNOTATE: a cell subset needs both a column and its values.")
    if subset_col not in adata.obs.columns:
        raise ValueError(
            f"ADATA_PREPANNOTATE: the subset column '{subset_col}' is not in obs. "
            f"Available: {short_list(adata.obs.columns)}"
        )
    wanted = list(dict.fromkeys(value.strip() for value in subset_values.split(",") if value.strip()))
    column = adata.obs[subset_col].astype(str)
    selected = column.isin(wanted).to_numpy()
    unmatched = [value for value in wanted if not (column == value).any()]
    if unmatched:
        warn(f"subset value(s) not found in '{subset_col}': {short_list(unmatched)}")
    if not selected.any():
        raise ValueError(
            f"ADATA_PREPANNOTATE: no cell has one of the values {wanted} in '{subset_col}'. "
            f"Values present: {short_list(sorted(column.unique()))}"
        )
    print_counts(f"Subset by '{subset_col}', cells per selected value", column[selected])
    print(f"Kept {int(selected.sum())} of {adata.n_obs} cells.")
    adata = adata[selected].copy()

# Raw counts, as a float32 CSR matrix

if counts_key not in adata.layers:
    raise ValueError(
        f"ADATA_PREPANNOTATE: the input has no layers['{counts_key}']. This stage needs the raw counts "
        f"there, as written by --stage qc_clustering. Available layers: {short_list(adata.layers.keys()) or 'none'}"
    )
counts = sp.csr_matrix(adata.layers[counts_key], dtype=np.float32)

values = counts.data
if not values.size or not np.any(values):
    raise ValueError(f"ADATA_PREPANNOTATE: layers['{counts_key}'] holds no counts (all zero).")
if not np.all(np.isfinite(values)):
    raise ValueError(f"ADATA_PREPANNOTATE: layers['{counts_key}'] has non-finite values (NaN or infinite).")
if np.any(values < 0):
    raise ValueError(f"ADATA_PREPANNOTATE: layers['{counts_key}'] has negative values; it must hold raw counts.")
is_integer = bool(np.all(values == np.round(values)))
max_value = float(values.max())
if not is_integer:
    if max_value <= max_lognorm_value:
        raise ValueError(
            f"ADATA_PREPANNOTATE: layers['{counts_key}'] looks log-normalised (non-integer values, maximum "
            f"{max_value:.2f}); it must hold raw counts."
        )
    # e.g. SoupX or decontX corrected counts
    warn(f"layers['{counts_key}'] has non-integer values (e.g. ambient RNA corrected counts); used as they are.")

# CellTypist cannot normalise, and scANVI cannot model, cells without counts
cell_totals = np.asarray(counts.sum(axis=1)).ravel()
n_empty = int((cell_totals == 0).sum())
if n_empty:
    raise ValueError(
        f"ADATA_PREPANNOTATE: {n_empty} cells have no counts in layers['{counts_key}']; remove them first."
    )

# obs: all columns (a scANVI model may need its covariates) except earlier predictions

obs = adata.obs.copy()
old_predictions = [col for col in obs.columns if str(col).startswith(prediction_prefixes)]
if old_predictions:
    print(f"Dropped {len(old_predictions)} earlier prediction columns from the working copy: {short_list(old_predictions)}")
    obs = obs.drop(columns=old_predictions)

# Batches (the scANVI query batches)

if batch_col and batch_col in obs.columns:
    batch_source = batch_col
elif fallback_batch_col in obs.columns:
    warn(f"the batch column '{batch_col}' is not in obs; using '{fallback_batch_col}' as the batches.")
    batch_source = fallback_batch_col
else:
    warn(
        f"neither the batch column '{batch_col}' nor '{fallback_batch_col}' is in obs; "
        f"the cells are mapped as a single batch '{single_batch}'."
    )
    batch_source = None

if batch_source is None:
    obs[batch_key] = pd.Categorical([single_batch] * len(obs))
else:
    if obs[batch_source].isna().any():
        warn(f"{int(obs[batch_source].isna().sum())} cells have no value in '{batch_source}'; they form the batch 'nan'.")
    obs[batch_key] = pd.Categorical(obs[batch_source].astype(str).to_numpy())

# Gene symbols

var = adata.var.copy()
if symbol_col not in ("index", ""):
    if symbol_col not in var.columns:
        raise ValueError(
            f"ADATA_PREPANNOTATE: the symbol column '{symbol_col}' is not in var. "
            f"Available: {short_list(var.columns)}"
        )
    symbols = var[symbol_col].astype(object)
    missing = symbols.isna() | (symbols.astype(str).str.strip() == "")
    if missing.any():
        warn(f"{int(missing.sum())} genes have no symbol in '{symbol_col}'; they keep their var_names.")
    if input_var_names_key not in var.columns:
        var[input_var_names_key] = var.index.astype(str)
    var.index = pd.Index(np.where(missing, var.index.astype(str), symbols.astype(str)))
var.index.name = None

adata = ad.AnnData(X=counts, obs=obs, var=var)
if not adata.var_names.is_unique:
    n_duplicated = int(adata.var_names.duplicated().sum())
    warn(f"{n_duplicated} duplicated gene names made unique by appending numbers.")
    adata.var_names_make_unique()
adata.layers[counts_key] = adata.X.copy()

adata.write_h5ad(f"{prefix}.h5ad")

print_counts(f"Query batches from '{batch_source or single_batch}', cells per batch", adata.obs[batch_key])
print(
    f"{adata.n_obs} cells, {adata.n_vars} genes, {adata.obs[batch_key].nunique()} batches; "
    f"counts are {'integers' if is_integer else 'non-integer'} (maximum {max_value:g})."
)

# Versions

versions = {
    "${task.process}": {
        "python": platform.python_version(),
        "anndata": ad.__version__,
        "numpy": np.__version__,
        "pandas": pd.__version__,
        "scipy": scipy.__version__,
    }
}

with open("versions.yml", "w") as f:
    yaml.dump(versions, f)
