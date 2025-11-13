#!/usr/bin/env python3

import os

os.environ["NUMBA_CACHE_DIR"] = "./tmp/numba"
os.environ["MPLCONFIGDIR"] = "./tmp/matplotlib"

import scanpy as sc
import matplotlib
import matplotlib.pyplot as plt
import platform
import base64
import json
import yaml


adata = sc.read_h5ad("${h5ad}")
prefix = "${prefix}"

## (currently) fixed parameters
dpi_plots = 120


# check for gene name, symbol etc.
for name in ['gene_name', 'gene_symbol']:
    if name in adata.var:
        adata.var.index = adata.var[name]
        break

# make var names unique
adata.var_names_make_unique()


## compute metrics and thresholds

# get mito, ribo, hemoglobin genes counts
adata.var["mt"] = adata.var_names.str.lower().str.startswith("mt-")
adata.var["ribo"] = adata.var_names.str.lower().str.startswith(("rps", "rpl"))
adata.var["hb"] = adata.var_names.str.lower().str.contains("^hb[^(p)]")

sc.pp.calculate_qc_metrics(
    adata, qc_vars=["mt", "ribo", "hb"], inplace=True, percent_top=None, log1p=True
)

# save anndata with QC metrics
adata.write_h5ad(f"{prefix}.h5ad")


## plots
sc.set_figure_params(dpi=dpi_plots, fontsize=18, dpi_save=dpi_plots)

# violin
sc.pl.violin(
    adata,
    ["n_genes_by_counts", "total_counts", "pct_counts_mt", "pct_counts_ribo", "pct_counts_hb"],
    multi_panel=True,
    stripplot=False
)
path_plt_violin = "${prefix}_qc_violin_plots.png"
plt.savefig(path_plt_violin)

# scatter
sc.pl.scatter(adata, x='total_counts', y='n_genes_by_counts', color="pct_counts_mt", show=False)
path_plt_scatter = "${prefix}_total_counts_vs_n_genes_by_counts.png"
plt.savefig(path_plt_scatter)


# MultiQC
with open("${prefix}_mqc.json", "w") as f_json:

    image_html = ""    

    with open(path_plt_violin, "rb") as f_plot:
        image_string = base64.b64encode(f_plot.read()).decode("utf-8")
        image_html += f'<figure><div class="mqc-custom-content-image"><img src="data:image/png;base64,{image_string}" /></div>'
        image_html += '<figcaption>The violin plots show the distribution of gene counts, total counts, percentage of mitochondrial (<i>mt</i>), ribosomal (<i>ribo</i>),'
        image_html += 'and hemoglobin gene expression (<i>hb</i>) per cell. Low gene and total counts may indicate poor quality cells '
        image_html += 'while high counts can be caused by artefacts. High mitochondrial, ribosmal and hemoglobin gene expression may indicate contamination, '
        image_html += 'degradation, and artefacts.</figcaption></figure>'

    with open(path_plt_scatter, "rb") as f_plot:
        image_string = base64.b64encode(f_plot.read()).decode("utf-8")
        image_html += f'<figure><div class="mqc-custom-content-image"><img src="data:image/png;base64,{image_string}" /></div>'
        image_html += '<figcaption>The scatter plot shows the relationship between total counts, gene counts,  mitochondrial(<i>mt</i>) gene expression.</figcaption></figure>'

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
