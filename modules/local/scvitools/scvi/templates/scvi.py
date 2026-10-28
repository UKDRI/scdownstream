#!/usr/bin/env python3

import os
import platform

os.environ["CUBLAS_WORKSPACE_CONFIG"] = ":4096:8"

import scvi
import anndata as ad
import pandas as pd
from scvi.model import SCVI
import torch
import yaml

from threadpoolctl import threadpool_limits
threadpool_limits(int("${task.cpus}"))

adata = ad.read_h5ad("${h5ad}")
prefix = "${prefix}"
mode = "${mode}"
reference_model_path = "reference_model"
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
# A model trained from scratch stops once the validation loss no longer improves
scratch_early_stopping = True
# scArches query mapping: no weight decay, so the frozen reference weights are not shrunk,
# and a fixed number of epochs without early stopping
query_weight_decay = 0.0
query_early_stopping = False

torch.set_float32_matmul_precision(matmul_precision)
torch.use_deterministic_algorithms(True)

scvi.settings.num_threads = int("${task.cpus}")
scvi.settings.seed = random_seed

# raw counts should be stored in 'counts' layer if not copy X
if "counts" not in adata.layers:
    adata.layers["counts"] = adata.X.copy()


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


x_input = None
added_layer = None

if mode == "query":
    if not os.path.exists(os.path.join(reference_model_path, "model.pt")):
        raise ValueError("mode 'query' needs a reference scVI model, but none was given.")

    registry = SCVI.load_registry(reference_model_path)
    if registry["model_name"] != "SCVI":
        raise ValueError(
            f"The reference model is a {registry['model_name']} model; SCVITOOLS_SCVI needs an SCVI model."
        )
    setup_args = registry["setup_args"]
    check_model_obs_columns(adata, setup_args, "reference scVI model")

    SCVI.prepare_query_anndata(adata, reference_model_path)

    # The model reads the layer it was set up with ("counts" for models trained by this pipeline).
    # Put the raw counts there; a model set up on X gets the counts in X while it runs.
    model_layer = setup_args.get("layer")
    if model_layer is None:
        x_input = adata.X
        adata.X = adata.layers["counts"].copy()
    elif model_layer not in adata.layers:
        adata.layers[model_layer] = adata.layers["counts"].copy()
        added_layer = model_layer

    model = SCVI.load_query_data(adata, reference_model_path)
    train_kwargs = {
        "max_epochs": query_max_epochs,
        "early_stopping": query_early_stopping,
        "plan_kwargs": {"weight_decay": query_weight_decay},
    }
elif mode == "scratch":
    categorical_covariates = (
        categorical_covariates.split(",") if categorical_covariates else None
    )
    continuous_covariates = (
        continuous_covariates.split(",") if continuous_covariates else None
    )

    SCVI.setup_anndata(
        adata,
        layer="counts",
        batch_key=batch_col,
        categorical_covariate_keys=categorical_covariates,
        continuous_covariate_keys=continuous_covariates,
    )

    model = SCVI(
        adata,
        n_hidden=int("${n_hidden}"),
        n_layers=int("${n_layers}"),
        n_latent=int("${n_latent}"),
        dispersion="${dispersion}",
        gene_likelihood="${gene_likelihood}",
    )
    train_kwargs = {
        "max_epochs": int(max_epochs) if max_epochs else None,
        "early_stopping": scratch_early_stopping,
    }
else:
    raise ValueError(f"Unknown mode '{mode}' for SCVITOOLS_SCVI; expected 'scratch' or 'query'.")

if use_gpu:
    model.to_device(0)

model.train(**train_kwargs)

latent = model.get_latent_representation()
adata.obsm["X_emb"] = latent
adata.obsm["X_scvi"] = latent

# Restore the input X / drop the helper layer, now that the model no longer reads them
if x_input is not None:
    adata.X = x_input
if added_layer is not None:
    del adata.layers[added_layer]

adata.uns.pop("_scvi_manager_uuid", None)
adata.uns.pop("_scvi_uuid", None)

adata.write_h5ad(f"{prefix}.h5ad")
model.save(f"{prefix}_model")

df = pd.DataFrame(adata.obsm["X_emb"], index=adata.obs_names)
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
