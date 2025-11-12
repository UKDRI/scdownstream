#!/usr/bin/env python3

import os

# These are needed to prevent errors during import of scanpy
# when using singularity/apptainer
os.environ["MPLCONFIGDIR"] = "./tmp/mpl"
os.environ["NUMBA_CACHE_DIR"] = "./tmp/numba"

import scanpy as sc
import platform
import yaml
import json
import base64
from threadpoolctl import threadpool_limits
threadpool_limits(int("${task.cpus}"))
sc.settings.n_jobs = int("${task.cpus}")

adata = sc.read_h5ad("${h5ad}")
prefix = "${prefix}"
batch_col = "${batch_col ?: ''}"

kwargs = {}
if batch_col and adata.obs[batch_col].nunique() > 1:
    kwargs["batch_key"] = batch_col

## (currently) fixed parameters
dpi_plots = 120

    
# run scrublet
sc.pp.scrublet(adata, **kwargs)

# plot doublet distribution
path_plt = "${prefix}_distributions.png"
fig = sc.pl.scrublet_score_distribution(adata, return_fig=True)
fig.savefig(path_plt)

# save anndata
adata.write_h5ad(f"{prefix}.h5ad")


adata.obs["predicted_doublet"] = adata.obs["predicted_doublet"].astype(bool)
df = adata.obs[["predicted_doublet"]]
df.columns = ["${prefix}"]
df.to_pickle("${prefix}.pkl")

# adata = adata[~adata.obs["predicted_doublet"]].copy()


# multiQC
with open("${prefix}_mqc.json", "w") as f_json:

    image_html = ""
    with open(path_plt, "rb") as f_plot:
        image_string = base64.b64encode(f_plot.read()).decode("utf-8")
        image_html += f'<figure><div class="mqc-custom-content-image"><img src="data:image/png;base64,{image_string}" /></div>'
        image_html += '<figcaption>You should see a bimodal distribution in the plots above. If there is no bimodal distribution the prediction of doublets is not realiable.'
        image_html += '</figcaption></figure>'
        
    custom_json = {
        "id": "${prefix}",
        "parent_id": "${section_name}".replace(" ", "_"),
        "parent_name": "${section_name}",
        "parent_description": "${description}",

        "section_name": "${meta.id}",
        "plot_type": "image",
        "data": image_html,
    }

    json.dump(custom_json, f_json)
        
# Versions

versions = {
    "${task.process}": {
        "python": platform.python_version(),
        "scanpy": sc.__version__
    }
}

with open("versions.yml", "w") as f:
    yaml.dump(versions, f)
