#!/usr/bin/env python3

import os
import platform

os.environ["CUBLAS_WORKSPACE_CONFIG"] = ":4096:8"

import scvi
import anndata as ad
import pandas as pd
from scvi.model import SCVI, SCANVI
import torch
import yaml

from threadpoolctl import threadpool_limits
threadpool_limits(int("${task.cpus}"))

adata = ad.read_h5ad("${h5ad}")
prefix = "${prefix}"
mode = "${mode}"
model_name = "${model_name}"
reference_model_path = "reference_model"
label_col = "${label_col}"
unlabeled_category = "${unlabeled_category}"
batch_col = "${batch_col}"
categorical_covariates = "${categorical_covariates}"
continuous_covariates = "${continuous_covariates}"
max_epochs = "${max_epochs ?: ''}"
query_max_epochs = int("${query_max_epochs}")
use_gpu = "${task.ext.use_gpu}" == "true"

# Fixed seed so that training is reproducible between runs
random_seed = 0
# "medium" allows bfloat16 matmuls on GPUs that support them (faster, slightly less precise)
matmul_precision = "medium"
# Training on labelled data (scratch / from_scvi) stops once the validation loss no longer improves
reference_early_stopping = True
# scArches query mapping: no weight decay, so the frozen reference weights are not shrunk,
# and a fixed number of epochs without early stopping
query_weight_decay = 0.0
query_early_stopping = False
# Name of the obs column that keeps the query's own labels, which are hidden from the model
label_input_col = "label_input"

torch.set_float32_matmul_precision(matmul_precision)
torch.use_deterministic_algorithms(True)

scvi.settings.num_threads = int("${task.cpus}")
scvi.settings.seed = random_seed

# raw counts should be stored in 'counts' layer if not copy X
if "counts" not in adata.layers:
    adata.layers["counts"] = adata.X.copy()

pred_col = f"scanvi:{model_name}"
conf_col = f"{pred_col}:conf"


def require_reference_model(mode):
    if not os.path.exists(os.path.join(reference_model_path, "model.pt")):
        raise ValueError(f"mode '{mode}' needs a reference model, but none was given.")


def load_checked_registry(model_cls, model_desc):
    registry = model_cls.load_registry(reference_model_path)
    if registry["model_name"] != model_cls.__name__:
        raise ValueError(
            f"The {model_desc} is a {registry['model_name']} model; mode '{mode}' needs a {model_cls.__name__} model."
        )
    return registry


def check_model_obs_columns(adata, setup_args, model_desc):
    # The query must have every obs column the reference model was set up with
    needed = [setup_args.get("batch_key"), setup_args.get("size_factor_key")]
    needed += setup_args.get("categorical_covariate_keys") or []
    needed += setup_args.get("continuous_covariate_keys") or []
    missing = [col for col in needed if col and col not in adata.obs.columns]
    if missing:
        raise ValueError(
            f"The {model_desc} needs the obs column(s) {missing}, which the input data does not have."
        )


def put_counts_in_model_layer(adata, model_layer):
    # The model reads the layer it was set up with ("counts" for models trained by this pipeline).
    # Put the raw counts there; a model set up on X gets the counts in X while it runs.
    # Returns what to undo once the model is done: (original X or None, added layer or None).
    if model_layer is None:
        x_input = adata.X
        adata.X = adata.layers["counts"].copy()
        return x_input, None
    if model_layer not in adata.layers:
        adata.layers[model_layer] = adata.layers["counts"].copy()
        return None, model_layer
    return None, None


def check_known_labels(adata, label_col, unlabeled_category):
    if label_col not in adata.obs.columns:
        raise ValueError(f"The label column '{label_col}' is not in the input data.")
    known_labels = set(adata.obs[label_col].dropna().unique())
    known_labels.discard(unlabeled_category)
    if len(known_labels) < 2:
        raise ValueError("Not enough labels to run scANVI.")


x_input, added_layer = None, None

if mode == "scratch":
    check_known_labels(adata, label_col, unlabeled_category)

    categorical_covariates = categorical_covariates.split(",") if categorical_covariates else None
    continuous_covariates = continuous_covariates.split(",") if continuous_covariates else None

    SCANVI.setup_anndata(
        adata,
        layer="counts",
        batch_key=batch_col,
        labels_key=label_col,
        unlabeled_category=unlabeled_category,
        categorical_covariate_keys=categorical_covariates,
        continuous_covariate_keys=continuous_covariates,
    )

    model = SCANVI(
        adata,
        n_hidden=int("${n_hidden}"),
        n_layers=int("${n_layers}"),
        n_latent=int("${n_latent}"),
        dispersion="${dispersion}",
        gene_likelihood="${gene_likelihood}",
    )
    train_kwargs = {
        "max_epochs": int(max_epochs) if max_epochs else None,
        "early_stopping": reference_early_stopping,
    }

elif mode == "from_scvi":
    require_reference_model(mode)
    check_known_labels(adata, label_col, unlabeled_category)
    registry = load_checked_registry(SCVI, "reference scVI model")
    check_model_obs_columns(adata, registry["setup_args"], "reference scVI model")

    SCVI.prepare_query_anndata(adata, reference_model_path)
    x_input, added_layer = put_counts_in_model_layer(adata, registry["setup_args"].get("layer"))

    scvi_model = SCVI.load(reference_model_path, adata)
    model = SCANVI.from_scvi_model(
        scvi_model,
        labels_key=label_col,
        unlabeled_category=unlabeled_category,
    )
    train_kwargs = {
        "max_epochs": int(max_epochs) if max_epochs else None,
        "early_stopping": reference_early_stopping,
    }

elif mode == "query":
    require_reference_model(mode)
    registry = load_checked_registry(SCANVI, "reference scANVI model")
    setup_args = registry["setup_args"]
    check_model_obs_columns(adata, setup_args, "reference scANVI model")

    # Hide the query's labels from the model: keep them as label_input, and mark every cell as
    # unlabelled with the model's own unlabelled category, so no unknown label category is introduced
    model_labels_key = setup_args["labels_key"]
    model_unlabeled_category = setup_args["unlabeled_category"]
    if model_labels_key in adata.obs.columns:
        adata.obs[label_input_col] = adata.obs[model_labels_key].copy()
    adata.obs[model_labels_key] = model_unlabeled_category

    SCANVI.prepare_query_anndata(adata, reference_model_path)
    x_input, added_layer = put_counts_in_model_layer(adata, setup_args.get("layer"))

    model = SCANVI.load_query_data(adata, reference_model_path)
    train_kwargs = {
        "max_epochs": query_max_epochs,
        "early_stopping": query_early_stopping,
        "plan_kwargs": {"weight_decay": query_weight_decay},
    }

else:
    raise ValueError(
        f"Unknown mode '{mode}' for SCVITOOLS_SCANVI; expected 'scratch', 'from_scvi' or 'query'."
    )

if use_gpu:
    model.to_device(0)

model.train(**train_kwargs)

latent = model.get_latent_representation()
adata.obsm["X_scanvi"] = latent
# kept for compatibility with consumers of the generic integration embedding
adata.obsm["X_emb"] = latent

# Label probabilities, cells x labels: the most likely label and its probability
probs = model.predict(soft=True)
adata.obs[pred_col] = probs.idxmax(axis=1).astype("category")
adata.obs[conf_col] = probs.max(axis=1)

# Restore the input X / drop the helper layer, now that the model no longer reads them
if x_input is not None:
    adata.X = x_input
if added_layer is not None:
    del adata.layers[added_layer]

adata.uns.pop("_scvi_manager_uuid", None)
adata.uns.pop("_scvi_uuid", None)

adata.write_h5ad(f"{prefix}.h5ad")
adata.obs[[pred_col, conf_col]].to_pickle(f"{prefix}.pkl")
model.save(f"{prefix}_model")

df = pd.DataFrame(adata.obsm["X_scanvi"], index=adata.obs_names)
df.to_pickle(f"X_{prefix}.pkl")

# Versions

versions = {
    "${task.process}": {
        "python": platform.python_version(),
        "scvi": scvi.__version__
    }
}

with open("versions.yml", "w") as f:
    yaml.dump(versions, f)
