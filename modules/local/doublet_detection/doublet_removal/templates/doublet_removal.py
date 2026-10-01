#!/opt/conda/bin/python

import os
import platform
import base64
import json

os.environ["MPLCONFIGDIR"] = "./tmp/matplotlib"

import anndata as ad
import pandas as pd
import matplotlib.pyplot as plt
import upsetplot
import matplotlib

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

adata = ad.read_h5ad("${h5ad}")
threshold = int("${threshold}")
prefix = "${prefix}"

# Smallest method intersection drawn in the upset plot, so it stays readable
min_subset_size = 10

def load(path: str) -> pd.DataFrame:
    if path.endswith(".pkl"):
        return pd.read_pickle(path)
    if path.endswith(".csv"):
        return pd.read_csv(path, index_col=0)

def method_name(path: str) -> str:
    # "<meta.id>_scrublet.pkl" -> "scrublet", otherwise the file stem
    stem = os.path.splitext(os.path.basename(path))[0]
    return stem.removeprefix("${meta.id}_") or stem

# One bool column per method, aligned to the cells of this object; cells without a prediction
# (not scored, or missing from a method's output) count as not doublets
predictions = pd.concat(
    [load(f).set_axis([method_name(f)], axis=1) for f in "${predictions}".split()], axis=1
).reindex(adata.obs_names)
for method, n_missing in predictions.isna().sum().items():
    if n_missing:
        print(f"WARNING: {n_missing} cells have no prediction from {method}; counted as not doublets")
predictions = predictions.astype("boolean").fillna(False).astype(bool)

n_methods = predictions.shape[1]
if threshold > n_methods:
    print(f"WARNING: threshold {threshold} is larger than the number of methods ({n_methods}); no cells are removed")

mask = predictions.sum(axis=1) >= threshold
n_cells = adata.n_obs
n_removed = int(mask.sum())
print(f"Removed {n_removed} of {n_cells} cells predicted as doublets by at least {threshold} of {n_methods} method(s)")

adata = adata[~mask.to_numpy()].copy()
adata.write_h5ad(f"{prefix}.h5ad")

# Versions

versions = {
    "${task.process}": {
        "python": platform.python_version(),
        "anndata": ad.__version__,
        "pandas": pd.__version__,
        "matplotlib": matplotlib.__version__,
        "upsetplot": upsetplot.__version__,
    }
}

with open("versions.yml", "w") as f:
    f.write(format_yaml_like(versions))

summary_html = (f"<p>{n_removed} of {n_cells} cells removed (predicted as doublets by at least "
                f"{threshold} of {n_methods} method(s)).</p>")

image_html = ""
plot_type = "html"

if n_methods > 1:
    # Plot

    contents = {column: predictions[column][predictions[column]].index.tolist()
                   for column in predictions.columns}

    plot_data = upsetplot.from_contents(contents)

    # cells per intersection of methods; upsetplot fails when no intersection reaches min_subset_size
    subset_sizes = plot_data.groupby(level=list(range(plot_data.index.nlevels))).size()

    if (subset_sizes >= min_subset_size).any():
        upsetplot.plot(plot_data,
                       sort_by="cardinality",
                       show_counts=True,
                       subset_size="count",
                       min_subset_size=min_subset_size)
        plot_path = f"{prefix}_predictions_mqc.png"
        plt.savefig(plot_path)

        with open(plot_path, "rb") as f_plot:
            image_string = base64.b64encode(f_plot.read()).decode("utf-8")
        image_html = f'<div class="mqc-custom-content-image"><img src="data:image/png;base64,{image_string}" /></div>'
        plot_type = "image"
    else:
        image_html = f"<p>No upset plot: every intersection of the methods has fewer than {min_subset_size} cells.</p>"

# MultiQC

with open("${prefix}_mqc.json", "w") as f_json:
    custom_json = {
        "id": "${prefix}",
        "parent_id": "doublet_predictions",
        "parent_name": "Doublet predictions",
        "parent_description": "Doublets removed per sample, and upset plots of the doublet prediction tools.",

        "section_name": "${meta.id}",
        "plot_type": plot_type,
        "data": summary_html + image_html,
    }

    json.dump(custom_json, f_json)
