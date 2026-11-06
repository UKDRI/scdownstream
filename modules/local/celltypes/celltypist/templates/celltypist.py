#!/usr/bin/env python3

import os
import platform
from pathlib import Path

os.environ["MPLCONFIGDIR"] = "./tmp/mpl"
os.environ["NUMBA_CACHE_DIR"] = "./tmp/numba"
os.environ["CELLTYPIST_FOLDER"] = "./tmp/celltypist"

import pandas as pd
import scanpy as sc
import celltypist
from celltypist import models as ct_models

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

models = "${models.join(' ')}".split()
# Local model files, staged by the process into this folder (stageAs: 'celltypist_models/*')
model_files = sorted(str(path) for path in Path("celltypist_models").glob("*") if path.is_file())

if not models and not model_files:
    raise ValueError("CELLTYPES_CELLTYPIST: no CellTypist model given (neither a model name nor a .pkl file).")

# (column name, model to load): built-in models are downloaded, local files are loaded as they are
model_specs = []
for model in models:
    model_file = f"{model}.pkl" if not model.endswith(".pkl") else model
    ct_models.download_models(model=model_file)
    model_specs.append((model_file[:-4], model_file))
for model_file in model_files:
    model_specs.append((Path(model_file).stem, model_file))

model_names = [name for name, _ in model_specs]
duplicated = sorted({name for name in model_names if model_names.count(name) > 1})
if duplicated:
    raise ValueError(f"CELLTYPES_CELLTYPIST: more than one model named {duplicated}; the prediction columns would clash.")

adata_celltypist = adata.copy()  # make a copy of our adata
sc.pp.normalize_per_cell(
    adata_celltypist, counts_per_cell_after=10**4
)  # normalize to 10,000 counts per cell
sc.pp.log1p(adata_celltypist)  # log-transform

symbol_col = "${symbol_col}"
if symbol_col != "index" and symbol_col:
    if symbol_col not in adata_celltypist.var.columns:
        raise ValueError(f"Symbol column {symbol_col} not found in adata.var.columns")
    adata_celltypist.var_names = adata_celltypist.var[symbol_col]

df_list = []

for model_name, model_file in model_specs:
    # A path with a '/' is opened as a file; a bare name is looked up in the downloaded models
    model_obj = ct_models.Model.load(model_file)

    predictions = celltypist.annotate(
        adata_celltypist, model=model_obj
    )
    predictions_adata = predictions.to_adata()

    df_celltypist = predictions_adata.obs.loc[
        adata.obs.index, ["predicted_labels", "conf_score"]
    ]

    df_celltypist.columns = [f"celltypist:{model_name}", f"celltypist:{model_name}:conf"]
    df_list.append(df_celltypist)

df_celltypist = pd.concat(df_list, axis=1)
df_celltypist.to_pickle("${prefix}.pkl")

adata.obs = pd.concat([adata.obs, df_celltypist], axis=1)
adata.write_h5ad(f"{prefix}.h5ad")

# Versions

versions = {
    "${task.process}": {
        "python": platform.python_version(),
        "pandas": pd.__version__,
        "scanpy": sc.__version__,
        "celltypist": celltypist.__version__
    }
}

with open("versions.yml", "w") as f:
    f.write(format_yaml_like(versions))
