#!/usr/bin/env python3

import os

os.environ["NUMBA_CACHE_DIR"] = "./tmp/numba"
os.environ["MPLCONFIGDIR"] = "./tmp/matplotlib"

import scanpy as sc
import seaborn as sns
import matplotlib
import matplotlib.pyplot as plt
import platform
import base64
import json
import yaml

adata = sc.read_h5ad("${h5ad}")
prefix = "${prefix}"

# check for gene name, symbol etc. -- TODO: add as a data preparation step before calling QC
for name in ['gene_name', 'gene_symbol']:
    if name in adata.var:
        adata.var.index = adata.var[name]
        break

# make var names unique
adata.var_names_make_unique()

# get mito, ribo, hemoglobin genes
adata.var["mt"] = adata.var_names.str.lower().str.startswith("mt-")
adata.var["ribo"] = adata.var_names.str.lower().str.startswith(("rps", "rpl"))
adata.var["hb"] = adata.var_names.str.lower().str.contains("^hb[^(p)]")

sc.pp.calculate_qc_metrics(
    adata, qc_vars=["mt", "ribo", "hb"], inplace=True, percent_top=None, log1p=True
)

# sc.pp.calculate_qc_metrics(adata, percent_top=None, log1p=False, inplace=True)

## save anndata with QC metrics
adata.write_h5ad(f"{prefix}.h5ad")


## plots
sc.pl.scatter(adata, x='total_counts', y='n_genes_by_counts', color="pct_counts_mt", show=False)
path_plt_counts = "${prefix}_total_counts_vs_n_genes_by_counts.png"
plt.savefig(path_plt_counts)

sc.pl.violin(adata, "pct_counts_mt")
path_plt_mt = "${prefix}_pct_counts_mt.png"
plt.savefig(path_plt_mt)

sc.pl.violin(adata, "pct_counts_ribo")
path_plt_ribo= "${prefix}_pct_counts_ribo.png"
plt.savefig(path_plt_ribo)

sc.pl.violin(adata, "pct_counts_hb")
path_plt_hb = "${prefix}_pct_counts_hb.png"
plt.savefig(path_plt_hb)

# MultiQC
with open("${prefix}_mqc.json", "w") as f_json:

    image_html = ""    
    for path in [path_plt_counts, path_plt_mt, path_plt_ribo, path_plt_hb]:
        with open(path, "rb") as f_plot:
            image_string = base64.b64encode(f_plot.read()).decode("utf-8")
            image_html += f'<div class="mqc-custom-content-image"><img src="data:image/png;base64,{image_string}" /></div>'
    
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
        "scanpy": sc.__version__,
        "matplotlib": matplotlib.__version__,
    }
}



with open("versions.yml", "w") as f:
    yaml.dump(versions, f)
