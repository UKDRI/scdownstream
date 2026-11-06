#!/usr/bin/env python3

import os
import platform

os.environ["MPLCONFIGDIR"] = "./tmp/mpl"
os.environ["NUMBA_CACHE_DIR"] = "./tmp/numba"
os.environ["CELLTYPIST_FOLDER"] = "./tmp/celltypist"

import numpy as np
import pandas as pd
import scanpy as sc
import sklearn
import celltypist

from threadpoolctl import threadpool_limits
threadpool_limits(int("${task.cpus}"))


def format_yaml_like(data: dict, indent: int = 0) -> str:
    """Formats a dictionary to a YAML-like string.

    Args:
        data (dict): The dictionary to format.
        indent (int): The current indentation level.

    Returns:
        str: A string formatted as YAML.
    """
    yaml_str = ""
    for key, value in data.items():
        spaces = "  " * indent
        if isinstance(value, dict):
            yaml_str += f"{spaces}{key}:\\n{format_yaml_like(value, indent + 1)}"
        else:
            yaml_str += f"{spaces}{key}: {value}\\n"
    return yaml_str


adata = sc.read_h5ad("${h5ad}")
prefix = "${prefix}"
feature_selection = "${feature_selection}" == "true"
use_sgd = "${use_sgd}" == "true"
mini_batch = "${mini_batch}" == "true"
batch_size = int("${batch_size}")
n_jobs = int("${task.cpus}")

# obs column, unlabelled category and counts layer written by ADATA_TRAINSPLIT
label_key = "label"
unknown_category = "unknown"
counts_key = "counts"
# CellTypist's input: log1p of counts normalised to this total per cell, exactly as
# CELLTYPES_CELLTYPIST prepares the data it annotates
counts_per_cell_after = 10**4
# Fixed seed for sklearn's solvers (the SGD pass of the feature selection shuffles the cells),
# so that the model is reproducible between runs
random_state = 0
# The data are normalised above; CellTypist's check of the first cell would only repeat that
check_expression = False
# CellTypist 1.6.3 builds LogisticRegression(multi_class='ovr') whenever SGD is off (also in the
# second pass of the feature selection); scikit-learn removed multi_class in this version
sklearn_without_multi_class = (1, 8)

sklearn_version = tuple(int(part) for part in sklearn.__version__.split(".")[:2] if part.isdigit())
if not use_sgd and sklearn_version >= sklearn_without_multi_class:
    raise ValueError(
        f"CELLTYPIST_TRAIN: CellTypist {celltypist.__version__}'s logistic-regression training needs "
        f"scikit-learn < 1.8, but this environment has {sklearn.__version__}. "
        "Train with SGD instead (--train_celltypist_use_sgd true)."
    )
if mini_batch and not use_sgd:
    raise ValueError(
        "CELLTYPIST_TRAIN: mini-batch training needs SGD (set --train_celltypist_use_sgd true "
        "with --train_celltypist_mini_batch true)."
    )

# CellTypist draws its mini-batches with numpy's global random generator, which it never seeds
np.random.seed(random_state)

# Unlabelled cells carry no label to learn
adata = adata[adata.obs[label_key].astype(str) != unknown_category].copy()
adata.obs[label_key] = pd.Categorical(adata.obs[label_key].astype(str))
n_labels = len(adata.obs[label_key].cat.categories)
if n_labels < 2:
    raise ValueError(f"CELLTYPIST_TRAIN: {n_labels} labels in '{label_key}'; at least 2 are needed.")

if counts_key in adata.layers:
    adata.X = adata.layers[counts_key]
for key in list(adata.layers.keys()):
    del adata.layers[key]

sc.pp.normalize_per_cell(adata, counts_per_cell_after=counts_per_cell_after)
sc.pp.log1p(adata)

print(f"Training on {adata.n_obs} cells, {adata.n_vars} genes and {n_labels} labels "
      f"(feature_selection={feature_selection}, use_SGD={use_sgd}, mini_batch={mini_batch}"
      f"{f', batch_size={batch_size}' if mini_batch else ''}).")

model = celltypist.train(
    adata,
    labels=label_key,
    check_expression=check_expression,
    n_jobs=n_jobs,
    use_SGD=use_sgd,
    mini_batch=mini_batch,
    batch_size=batch_size,
    feature_selection=feature_selection,
    # passed on to sklearn's LogisticRegression / SGDClassifier
    random_state=random_state,
)

model.write(f"{prefix}.pkl")

# Versions

versions = {
    "${task.process}": {
        "python": platform.python_version(),
        "pandas": pd.__version__,
        "scanpy": sc.__version__,
        "celltypist": celltypist.__version__,
        "scikit-learn": sklearn.__version__,
    }
}

with open("versions.yml", "w") as f:
    f.write(format_yaml_like(versions))
