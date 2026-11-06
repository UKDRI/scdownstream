#!/usr/bin/env python3

import platform

import anndata as ad
import numpy as np
import pandas as pd
import sklearn
from sklearn.metrics import (
    accuracy_score,
    balanced_accuracy_score,
    confusion_matrix,
    f1_score,
    precision_recall_fscore_support,
)
import yaml

prefix = "${prefix}"
pred_col = "${pred_col}"
true_col = "${true_col}"

# Only obs is needed; backed mode leaves X on disk
obs = ad.read_h5ad("${h5ad}", backed="r").obs

# True labels with this value are unlabelled cells and are not evaluated (written by ADATA_TRAINSPLIT)
unknown_category = "unknown"
# A cell without a prediction counts as wrong, under this name
missing_prediction = "no prediction"
# Metrics that are 0/0 (e.g. the precision of a label that is never predicted) are set to 0 by
# sklearn and replaced by NaN in the per-label table
zero_division = 0

for col in (pred_col, true_col):
    if col not in obs.columns:
        raise ValueError(f"MODEL_EVALUATE: the column '{col}' is not in obs. Available: {', '.join(map(str, obs.columns))}")

has_truth = obs[true_col].notna().to_numpy() & (obs[true_col].astype(str) != unknown_category).to_numpy()
if not has_truth.any():
    raise ValueError(f"MODEL_EVALUATE: no cells with a true label in '{true_col}'.")

y_true = obs.loc[has_truth, true_col].astype(str).to_numpy()
y_pred = obs.loc[has_truth, pred_col].astype(object).where(obs.loc[has_truth, pred_col].notna(), missing_prediction)
y_pred = y_pred.astype(str).to_numpy()

# Rows of the confusion matrix and the labels the averages run over: the true labels.
# Predicted labels that are not in the truth are extra columns, and count as errors.
true_labels = sorted(np.unique(y_true))
extra_labels = sorted(set(np.unique(y_pred)) - set(true_labels))
all_labels = true_labels + extra_labels

# Overall metrics; macro and weighted F1 average over the true labels
metrics = pd.DataFrame(
    [
        ("n_cells", len(y_true)),
        ("accuracy", accuracy_score(y_true, y_pred)),
        ("balanced_accuracy", balanced_accuracy_score(y_true, y_pred)),
        ("macro_f1", f1_score(y_true, y_pred, labels=true_labels, average="macro", zero_division=zero_division)),
        ("weighted_f1", f1_score(y_true, y_pred, labels=true_labels, average="weighted", zero_division=zero_division)),
    ],
    columns=["metric", "value"],
)
metrics.to_csv(f"{prefix}_metrics.tsv", sep="\\t", index=False)

# Per label
precision, recall, f1, support = precision_recall_fscore_support(
    y_true, y_pred, labels=all_labels, average=None, zero_division=zero_division
)
n_predicted = pd.Series(y_pred).value_counts().reindex(all_labels, fill_value=0).to_numpy()
per_label = pd.DataFrame({
    "label": all_labels,
    "precision": precision,
    "recall": recall,
    "f1": f1,
    "support": support,
    "n_predicted": n_predicted,
})
# Undefined values: precision of a label never predicted, recall and F1 of a label not in the truth
per_label.loc[per_label["n_predicted"] == 0, "precision"] = np.nan
per_label.loc[per_label["support"] == 0, ["recall", "f1"]] = np.nan
per_label.to_csv(f"{prefix}_per_label.tsv", sep="\\t", index=False)

# Confusion matrix: counts, rows true labels, columns predicted labels
confusion = confusion_matrix(y_true, y_pred, labels=all_labels)[: len(true_labels)]
confusion = pd.DataFrame(confusion, index=pd.Index(true_labels, name="true_label"), columns=all_labels)
confusion.to_csv(f"{prefix}_confusion.tsv", sep="\\t")

print(metrics.to_string(index=False))

# Versions

versions = {
    "${task.process}": {
        "python": platform.python_version(),
        "anndata": ad.__version__,
        "numpy": np.__version__,
        "pandas": pd.__version__,
        "scikit-learn": sklearn.__version__,
    }
}

with open("versions.yml", "w") as f:
    yaml.dump(versions, f)
