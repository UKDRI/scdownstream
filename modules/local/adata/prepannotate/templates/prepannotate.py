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
prediction_cols = "${prediction_cols}"

# Names the annotation modules read (CELLTYPES_CELLTYPIST, SCVITOOLS_SCANVI with the batch key of the
# scANVI models); they define the structure of the output, not a setting
counts_key = "counts"
# Batch key of the scANVI models: SCVITOOLS_SCANVI reads the query batches from obs['batch']
batch_key = "batch"

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
    warn(f"the input has no layers['{counts_key}']; it will be created from X.")
    adata.layers[counts_key] = adata.X.copy()
    
if not sp.issparse(adata.layers[counts_key]):
    warn(f"the input layers['{counts_key}'] is not sparse; it will be converted to a sparse matrix.")
    adata.layers[counts_key] = sp.csr_matrix(adata.layers[counts_key], dtype=np.float32)

if not adata.layers[counts_key].size or not np.any(adata.layers[counts_key].data):
    raise ValueError(f"ADATA_PREPANNOTATE: layers['{counts_key}'] holds no counts (all zero).")
if not np.all(np.isfinite(adata.layers[counts_key].data)):
    raise ValueError(f"ADATA_PREPANNOTATE: layers['{counts_key}'] has non-finite values (NaN or infinite).")
if np.any(adata.layers[counts_key].data < 0):
    raise ValueError(f"ADATA_PREPANNOTATE: layers['{counts_key}'] has negative values; it must hold raw counts.")

# CellTypist cannot normalise, and scANVI cannot model, cells without counts
cell_totals = np.asarray(adata.layers[counts_key].sum(axis=1)).ravel()
n_empty = int((cell_totals == 0).sum())
if n_empty:
    raise ValueError(
        f"ADATA_PREPANNOTATE: {n_empty} cells have no counts in layers['{counts_key}']; remove them first."
    )

# Earlier predictions of the models this run predicts with are dropped, so the annotation modules
# only return the new ones (they would otherwise return the old column next to / instead of the new
# one); predictions of other models stay untouched
prediction_col_list = list(dict.fromkeys(col.strip() for col in prediction_cols.split(",") if col.strip()))
old_predictions = [
    name for col in prediction_col_list for name in (col, f"{col}:conf") if name in adata.obs.columns
]
# SCANVI_ANNOTATE latent key, X_scanvi_<model>; it does not exist for the CellTypist columns
old_latents = [key for key in ("X_" + col.replace(":", "_", 1) for col in prediction_col_list) if key in adata.obsm]
if old_predictions or old_latents:
    adata.obs = adata.obs.drop(columns=old_predictions)
    for key in old_latents:
        del adata.obsm[key]
    print(
        f"Dropped {len(old_predictions) + len(old_latents)} earlier prediction columns from the working copy: "
        f"{short_list(old_predictions + old_latents)}"
    )

# if no batches add all as batch 1
if batch_col and batch_col in adata.obs.columns:
    adata.obs[batch_key] = pd.Categorical(adata.obs[batch_col].astype(str).to_numpy())
else:
    warn(f"did not find a batch column '{batch_col}' in obs; all cells are mapped as a single batch '1'.")
    adata.obs[batch_key] = pd.Categorical(["1"] * adata.n_obs)

# Gene symbols
if symbol_col not in ("index", ""):
    if symbol_col not in adata.var.columns:
        raise ValueError(
            f"ADATA_PREPANNOTATE: the symbol column '{symbol_col}' is not in var. "
            f"Available: {short_list(adata.var.columns)}"
        )
    symbols = adata.var[symbol_col].astype(object)
    missing = symbols.isna() | (symbols.astype(str).str.strip() == "")
    if missing.any():
        warn(f"{int(missing.sum())} genes have no symbol in '{symbol_col}'; they keep their var_names.")
    if input_var_names_key not in adata.var.columns:
        adata.var[input_var_names_key] = adata.var.index.astype(str)
    adata.var.index = pd.Index(np.where(missing, adata.var.index.astype(str), symbols.astype(str)))

adata.var.index.name = None

if not adata.var_names.is_unique:
    n_duplicated = int(adata.var_names.duplicated().sum())
    warn(f"{n_duplicated} duplicated gene names made unique by appending numbers.")
    adata.var_names_make_unique()

# X is kept as in the input; the annotation modules take the raw counts from layers['counts']
adata.write_h5ad(f"{prefix}.h5ad")

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
