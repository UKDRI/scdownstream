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
import sklearn
from sklearn.model_selection import StratifiedGroupKFold
import yaml

from threadpoolctl import threadpool_limits
threadpool_limits(int("${task.cpus}"))

adata = ad.read_h5ad("${h5ad}")
prefix = "${prefix}"
label_col = "${label_col}"
unknown_label = "${unknown_label}"
batch_col = "${batch_col}"
counts_layer = "${counts_layer}"
symbol_col = "${symbol_col}"
min_cells_per_label = int("${min_cells_per_label}")
test_fraction = float("${test_fraction}")
split_by = "${split_by}"
evaluate = "${evaluate}" == "true"
# 0: no cap
max_cells_per_label = int("${max_cells_per_label}")

# Fixed seed so that the train/test split is reproducible between runs
random_seed = 0
# Fixed seed so that the downsampling to max_cells_per_label is reproducible between runs
downsample_seed = 0

# Names the training modules downstream read (CELLTYPIST_TRAIN, SCVITOOLS_*, MODEL_EVALUATE);
# they define the structure of the output, not a setting
label_key = "label"
true_label_key = "label_true"
batch_key = "batch"
counts_key = "counts"
unknown_category = "unknown"
# Batch of every cell when the input has no batch column
single_batch = "batch_1"
# var column that keeps the input var_names when the gene symbols come from another column
input_var_names_key = "var_names_input"


def warn(message):
    print(f"WARNING: ADATA_TRAINSPLIT: {message}", file=sys.stderr)


def short_list(values, n=20):
    values = [str(v) for v in values]
    return ", ".join(values[:n]) + (f", ... ({len(values)} in total)" if len(values) > n else "")


# Checks

if not label_col:
    raise ValueError("ADATA_TRAINSPLIT: no label column given (set --train_label_col).")
if label_col not in adata.obs.columns:
    raise ValueError(
        f"ADATA_TRAINSPLIT: the label column '{label_col}' is not in obs. "
        f"Available: {short_list(adata.obs.columns)}"
    )
if evaluate and not 0 < test_fraction < 1:
    raise ValueError(f"ADATA_TRAINSPLIT: the test fraction must be between 0 and 1, got {test_fraction}.")
if evaluate and split_by and split_by not in adata.obs.columns:
    raise ValueError(
        f"ADATA_TRAINSPLIT: the split column '{split_by}' (--train_test_split_by) is not in obs. "
        f"Available: {short_list(adata.obs.columns)}"
    )

# Raw counts, as a float32 CSR matrix

if counts_layer in ("X", ""):
    counts = adata.X
elif counts_layer in adata.layers:
    counts = adata.layers[counts_layer]
else:
    raise ValueError(
        f"ADATA_TRAINSPLIT: the counts layer '{counts_layer}' is not in the input. "
        f"Use 'X' or one of the layers: {short_list(adata.layers.keys())}"
    )
if counts is None:
    raise ValueError(f"ADATA_TRAINSPLIT: the counts matrix '{counts_layer}' is empty.")
counts = sp.csr_matrix(counts, dtype=np.float32)

if counts.data.size:
    if np.any(counts.data < 0):
        warn(f"'{counts_layer}' has negative values; the models expect raw counts (set --train_counts_layer).")
    elif np.any(counts.data != np.round(counts.data)):
        warn(f"'{counts_layer}' has non-integer values; the models expect raw counts (set --train_counts_layer).")

# Gene symbols

var = adata.var.copy()
if symbol_col not in ("index", ""):
    if symbol_col not in var.columns:
        raise ValueError(
            f"ADATA_TRAINSPLIT: the symbol column '{symbol_col}' is not in var. "
            f"Available: {short_list(var.columns)}"
        )
    symbols = var[symbol_col]
    n_missing = int(symbols.isna().sum())
    if n_missing:
        warn(f"{n_missing} genes have no symbol in '{symbol_col}'.")
    if input_var_names_key not in var.columns:
        var[input_var_names_key] = var.index.astype(str)
    var.index = pd.Index(symbols.astype(str).to_numpy())
var.index.name = None

# Labels, exactly as given; only the unknown label and missing values become unknown_category

obs = adata.obs.copy()
raw_labels = obs[label_col]
labels = raw_labels.astype(str)
unlabelled = (
    raw_labels.isna().to_numpy()
    | (labels.str.strip() == "").to_numpy()
    | (labels == unknown_label).to_numpy()
)
if unknown_label != unknown_category and (labels[~unlabelled] == unknown_category).any():
    warn(f"the label '{unknown_category}' is also treated as unlabelled.")
labels = labels.where(~unlabelled, unknown_category)
obs[label_key] = labels.to_numpy(dtype=object)

# Batches

if batch_col and batch_col in obs.columns:
    if obs[batch_col].isna().any():
        warn(f"{int(obs[batch_col].isna().sum())} cells have no value in '{batch_col}'; they form the batch 'nan'.")
    obs[batch_key] = obs[batch_col].astype(str).to_numpy(dtype=object)
else:
    warn(f"the batch column '{batch_col}' is not in obs; all cells are one batch '{single_batch}'.")
    obs[batch_key] = single_batch

adata = ad.AnnData(X=counts, obs=obs, var=var, obsm=dict(adata.obsm))
adata.var_names_make_unique()
if not adata.obs_names.is_unique:
    warn("cell names are not unique; made unique by appending numbers.")
    adata.obs_names_make_unique()

# Cells without any counts cannot be normalised (CellTypist) or modelled (scVI)
has_counts = np.asarray(adata.X.sum(axis=1)).ravel() > 0
if not has_counts.all():
    warn(f"removed {int((~has_counts).sum())} cells with no counts.")
    adata = adata[has_counts].copy()

adata.layers[counts_key] = adata.X.copy()

# Labels with too few cells are dropped everywhere, so the final model never knows a class
# the evaluation did not test

is_known = (adata.obs[label_key] != unknown_category).to_numpy()
label_sizes = adata.obs.loc[is_known, label_key].value_counts()
dropped_labels = sorted(label_sizes.index[label_sizes < min_cells_per_label])
kept_labels = sorted(label_sizes.index[label_sizes >= min_cells_per_label])

if dropped_labels:
    print(
        f"Dropped {len(dropped_labels)} labels with fewer than {min_cells_per_label} cells: "
        + "; ".join(f"{label} ({label_sizes[label]})" for label in dropped_labels)
    )
if len(kept_labels) < 2:
    raise ValueError(
        f"ADATA_TRAINSPLIT: {len(kept_labels)} labels with at least {min_cells_per_label} cells remain "
        f"in '{label_col}'; at least 2 are needed (lower --train_min_cells_per_label?)."
    )

adata = adata[~adata.obs[label_key].isin(dropped_labels).to_numpy()].copy()
n_unknown = int((adata.obs[label_key] == unknown_category).sum())

# Split

is_test = np.zeros(adata.n_obs, dtype=bool)
if evaluate:
    label_values = adata.obs[label_key].to_numpy(dtype=object)
    # Unlabelled cells are never tested; they stay in the training split
    labelled_pos = np.flatnonzero(label_values != unknown_category)
    labelled_y = label_values[labelled_pos]

    if split_by:
        groups = adata.obs[split_by].astype(str).to_numpy()[labelled_pos]
        n_groups = len(np.unique(groups))
        n_splits = int(round(1 / test_fraction))
        if n_splits < 2:
            raise ValueError(
                f"ADATA_TRAINSPLIT: a group split needs a test fraction below 2/3, got {test_fraction}."
            )
        if n_groups < n_splits:
            raise ValueError(
                f"ADATA_TRAINSPLIT: '{split_by}' has {n_groups} groups among the labelled cells; holding out "
                f"a fraction of {test_fraction} needs at least {n_splits}."
            )
        splitter = StratifiedGroupKFold(n_splits=n_splits, shuffle=True, random_state=random_seed)
        _, test_idx = next(splitter.split(np.zeros((len(labelled_pos), 1)), labelled_y, groups))
        is_test[labelled_pos[test_idx]] = True
        print(
            f"Group split by '{split_by}': {len(np.unique(groups[test_idx]))} of {n_groups} groups held out "
            f"(first of {n_splits} folds)."
        )
    else:
        rng = np.random.default_rng(random_seed)
        for label in kept_labels:
            pos = labelled_pos[labelled_y == label]
            # at least one cell stays in training
            n_test = min(int(np.floor(test_fraction * len(pos) + 0.5)), len(pos) - 1)
            if n_test > 0:
                is_test[rng.choice(pos, size=n_test, replace=False)] = True

# Downsampling: at most max_cells_per_label cells of each label (unlabelled cells count as one
# label) in the training and the all output, drawn independently; the test set is never capped


def cap_per_label(label_values, max_cells, seed):
    """Boolean mask keeping at most max_cells random cells of every label (all cells if max_cells <= 0)."""
    keep = np.ones(len(label_values), dtype=bool)
    if max_cells <= 0:
        return keep
    rng = np.random.default_rng(seed)
    for label in sorted(set(label_values)):
        pos = np.flatnonzero(label_values == label)
        if len(pos) > max_cells:
            keep[pos] = False
            keep[rng.choice(pos, size=max_cells, replace=False)] = True
    return keep


all_label_values = adata.obs[label_key].to_numpy(dtype=object)
keep_all = cap_per_label(all_label_values, max_cells_per_label, downsample_seed)
train_pos = np.flatnonzero(~is_test)
keep_train = np.zeros(adata.n_obs, dtype=bool)
keep_train[train_pos[cap_per_label(all_label_values[train_pos], max_cells_per_label, downsample_seed)]] = True

if max_cells_per_label > 0:
    print(
        f"Capped at {max_cells_per_label} cells per label: all {adata.n_obs} -> {int(keep_all.sum())} cells"
        + (f", train {len(train_pos)} -> {int(keep_train.sum())} cells" if evaluate else "")
        + "; the test set is not capped."
    )

# Label counts

n_train = adata.obs.loc[~is_test, label_key].value_counts()
n_test = adata.obs.loc[is_test, label_key].value_counts()
n_all_used = adata.obs.loc[keep_all, label_key].value_counts()
# Without evaluation there is no training split: the models are trained on the all output
n_train_used = adata.obs.loc[keep_train, label_key].value_counts() if evaluate else n_all_used

rows = []
for label in kept_labels + dropped_labels:
    dropped = label in dropped_labels
    rows.append({
        "label": label,
        "n_total": int(label_sizes[label]),
        "n_train": 0 if dropped else int(n_train.get(label, 0)),
        "n_test": 0 if dropped else int(n_test.get(label, 0)),
        "dropped": dropped,
        "n_train_used": 0 if dropped else int(n_train_used.get(label, 0)),
        "n_all_used": 0 if dropped else int(n_all_used.get(label, 0)),
    })
label_counts_columns = ["label", "n_total", "n_train", "n_test", "dropped", "n_train_used", "n_all_used"]
label_counts = pd.DataFrame(rows, columns=label_counts_columns)
label_counts = label_counts.sort_values(["dropped", "n_total", "label"], ascending=[True, False, True])
if n_unknown:
    label_counts = pd.concat([label_counts, pd.DataFrame([{
        "label": unknown_category, "n_total": n_unknown, "n_train": n_unknown, "n_test": 0, "dropped": False,
        "n_train_used": int(n_train_used.get(unknown_category, 0)),
        "n_all_used": int(n_all_used.get(unknown_category, 0)),
    }], columns=label_counts_columns)], ignore_index=True)
label_counts.to_csv(f"{prefix}_label_counts.tsv", sep="\\t", index=False)

if evaluate:
    kept_counts = label_counts[~label_counts["dropped"] & (label_counts["label"] != unknown_category)]
    no_test = kept_counts.loc[kept_counts["n_test"] == 0, "label"].tolist()
    no_train = kept_counts.loc[kept_counts["n_train"] == 0, "label"].tolist()
    if no_test:
        warn(f"{len(no_test)} labels have no test cells and are not evaluated: {short_list(no_test)}")
    if no_train:
        warn(f"{len(no_train)} labels have no training cells, so the evaluated models cannot predict them: "
             f"{short_list(no_train)}")

# Outputs


def with_label_categories(adata_out, key):
    adata_out.obs[key] = pd.Categorical(adata_out.obs[key].astype(str)).remove_unused_categories()
    return adata_out


def subset(mask):
    # no copy when nothing is removed
    return adata if mask.all() else adata[mask].copy()


with_label_categories(subset(keep_all), label_key).write_h5ad(f"{prefix}_all.h5ad")

if evaluate:
    with_label_categories(adata[keep_train].copy(), label_key).write_h5ad(f"{prefix}_train.h5ad")

    adata_test = adata[is_test].copy()
    adata_test.obs[true_label_key] = adata_test.obs[label_key].astype(str)
    adata_test.obs[label_key] = unknown_category
    adata_test = with_label_categories(with_label_categories(adata_test, true_label_key), label_key)
    adata_test.write_h5ad(f"{prefix}_test.h5ad")

print(
    f"{adata.n_obs} cells, {adata.n_vars} genes, {len(kept_labels)} labels kept, {len(dropped_labels)} dropped, "
    f"{n_unknown} unlabelled cells; train {int((~is_test).sum())}, test {int(is_test.sum())} cells "
    f"(written: all {int(keep_all.sum())}, train {int(keep_train.sum()) if evaluate else 0}, "
    f"test {int(is_test.sum()) if evaluate else 0})."
)

# Versions

versions = {
    "${task.process}": {
        "python": platform.python_version(),
        "anndata": ad.__version__,
        "numpy": np.__version__,
        "pandas": pd.__version__,
        "scipy": scipy.__version__,
        "scikit-learn": sklearn.__version__,
    }
}

with open("versions.yml", "w") as f:
    yaml.dump(versions, f)
