#!/usr/bin/env python3

import os

os.environ["NUMBA_CACHE_DIR"] = "./tmp/numba"

import anndata as ad
import scipy
import numpy as np
from scipy.sparse import csr_matrix
import platform
import yaml

def aggregate_duplicate_var(adata: ad.AnnData, method: str = "mean") -> ad.AnnData:
    """Collapse columns that share a var name into one column per unique name.

    Columns come out in first-occurrence order and keep the var row of the first
    occurrence. "sum" adds the duplicate columns, "mean" averages them (implicit zeros
    included) and "max" takes their element-wise maximum. X stays a CSR matrix of the
    same dtype and is never densified.
    """
    is_first = ~adata.var_names.duplicated(keep="first")
    if is_first.all():
        return adata

    X = csr_matrix(adata.X)
    dtype = X.dtype
    n_obs, n_var = X.shape
    first_idx = np.flatnonzero(is_first)
    n_unique = len(first_idx)

    # Column of each var in the deduplicated matrix. get_indexer matches NaN to NaN,
    # consistent with duplicated() above.
    codes = adata.var_names[first_idx].get_indexer(adata.var_names)
    assert (codes >= 0).all(), "Could not map every var name to its first occurrence"
    group_size = np.bincount(codes, minlength=n_unique)

    if method in ("sum", "mean"):
        # Indicator matrix: entry (j, k) is 1 if var j belongs to unique var k
        indicator = csr_matrix(
            (np.ones(n_var, dtype=dtype), (np.arange(n_var), codes)),
            shape=(n_var, n_unique),
        )
        X_new = (X @ indicator).tocsr()
        if method == "mean":
            X_new.data = (X_new.data / group_size[X_new.indices]).astype(dtype)
    elif method == "max":
        # Columns of unique vars are already correct in the first-occurrence slice;
        # only the duplicated groups need their maximum computed.
        is_dup_group = group_size > 1
        dup_cols = np.flatnonzero(is_dup_group[codes])
        dup = X[:, dup_cols].tocoo()
        rows, cols, vals = dup.row, codes[dup_cols][dup.col], dup.data

        # Sort stored entries by (row, group) and reduce each run to its maximum
        order = np.lexsort((cols, rows))
        rows, cols, vals = rows[order], cols[order], vals[order]
        if vals.size:
            starts = np.flatnonzero(np.r_[True, (rows[1:] != rows[:-1]) | (cols[1:] != cols[:-1])])
            group_max = np.maximum.reduceat(vals, starts)
            n_stored = np.diff(np.r_[starts, vals.size])
            # A run with fewer stored entries than its group has columns also covers
            # implicit zeros, which take part in the maximum
            has_implicit_zero = n_stored < group_size[cols[starts]]
            group_max = np.where(has_implicit_zero, np.maximum(group_max, 0), group_max)
            rows, cols = rows[starts], cols[starts]
        else:
            group_max = vals
        dup_max = csr_matrix((group_max, (rows, cols)), shape=(n_obs, n_unique), dtype=dtype)

        # Clear the duplicated groups' columns in the first-occurrence slice, then add the maxima
        X_first = X[:, first_idx]
        X_first.data[is_dup_group[X_first.indices]] = 0
        X_new = (X_first + dup_max).tocsr()
    else:
        raise ValueError(f"Invalid aggregation method: {method}")

    X_new = X_new.astype(dtype, copy=False)
    X_new.eliminate_zeros()

    adata_dedup = adata[:, is_first].copy()
    adata_dedup.X = X_new
    return adata_dedup

def to_Florent_case(s: str):
    corrected = s.lower().strip()

    if corrected in ["na", "nan", "null", "unknown"]:
        return "unknown"

    corrected = s \
        .replace(" ", "_") \
        .replace("-", "_")

    corrected = "".join([c if c.isalnum() or c == "_" else "" for c in corrected])

    # Make sure there is never more than one underscore
    corrected = corrected.replace("__", "_")

    if corrected.endswith("s"):
        corrected = corrected[:-1]

    corrected = corrected.strip(" _")

    if not corrected:
        return "unknown"

    return corrected.capitalize()

adata = ad.read_h5ad("$h5ad")

counts_layer = "${counts_layer}"
if counts_layer != "X":
    adata.X = adata.layers[counts_layer]

# Remove all obsm, varm, uns and layers
adata.obsm = {}
adata.varm = {}
adata.uns = {}
adata.layers = {}

# Convert to float32 CSR matrix
adata.X = csr_matrix(adata.X.astype(np.float32))

# Unify gene symbols
symbol_col = "${symbol_col}"

if symbol_col != "index":
    assert symbol_col in adata.var.columns, f"Symbol column {symbol_col} not found in var table"
    adata.var["original_index"] = adata.var.index
    adata.var.index = adata.var[symbol_col]

if "${aggregate_isoforms}" == "true":
    # Remove all numeric suffixes following a dot, keep non-numeric suffixes
    adata.var_names = adata.var_names.str.replace(r'\\.\\d+', '', regex=True)

# Deal with duplicate genes
method = "${duplicate_var_resolution}"
if method in ["mean", "sum", "max"]:
    adata = aggregate_duplicate_var(adata, method=method)
elif method == "make_unique":
    adata.var_names_make_unique()
else:
    raise ValueError(f"Invalid aggregation method: {method}")

# Prevent duplicate cells
adata.obs_names_make_unique()
adata.obs_names = "${meta.id}_" + adata.obs_names

# Unify batches
batch_col = "${batch_col}"
if batch_col not in adata.obs:
    adata.obs[batch_col] = "${meta.id}"

if batch_col != "batch":
    if "batch" in adata.obs:
        raise ValueError("The batch column already exists.")
    adata.obs["batch"] = adata.obs[batch_col]
    del adata.obs[batch_col]
adata.obs["batch"] = adata.obs["batch"].astype(str).astype("category")

# Unify labels
label_col = "${label_col}"
unknown_label = "${unknown_label}"

if label_col:
    if label_col not in adata.obs:
        raise ValueError("The specified label column does not exist in the dataset. Existing columns: " + ", ".join(adata.obs.columns))

    if label_col != "label":
        if "label" in adata.obs:
            raise ValueError("The label column already exists.")
        adata.obs["label"] = adata.obs[label_col]
        del adata.obs[label_col]

    if unknown_label != "unknown":
        if "unknown" in adata.obs["label"]:
            raise ValueError("The label column already contains 'unknown' values.")
        adata.obs["label"].replace({unknown_label: "unknown"}, inplace=True)

    # Replace all NaN values with "unknown"
    adata.obs["label"] = adata.obs["label"].astype(str)
    adata.obs["label"] = adata.obs["label"].fillna("unknown")
    adata.obs["label"] = adata.obs["label"].map(to_Florent_case)
    adata.obs["label"] = adata.obs["label"].astype("category")
else:
    if "label" in adata.obs:
        raise ValueError("The label column already exists.")
    adata.obs["label"] = "unknown"
adata.obs["label"] = adata.obs["label"].astype("category")

# Add "sample" column, keeping a differing pre-existing one as "sample_original"
if "sample" in adata.obs and not (adata.obs["sample"].astype(str) == "${meta.id}").all():
    if "sample_original" in adata.obs:
        raise ValueError("obs['sample'] holds values other than '${meta.id}', but obs['sample_original'] already exists, so it cannot be kept.")
    adata.obs["sample_original"] = adata.obs["sample"]
adata.obs["sample"] = "${meta.id}"
adata.obs["sample"] = adata.obs["sample"].astype("category")

# Add sample to batch column, to avoid overlap with other samples
adata.obs["batch"] = adata.obs["batch"].astype(str) + "_" + adata.obs["sample"].astype(str)
adata.obs["batch"] = adata.obs["batch"].astype("category")

adata.write_h5ad("${prefix}.h5ad")

# Versions

versions = {
    "${task.process}": {
        "python": platform.python_version(),
        "anndata": ad.__version__,
        "scipy": scipy.__version__,
        "numpy": np.__version__
    }
}

with open("versions.yml", "w") as f:
    yaml.dump(versions, f)
