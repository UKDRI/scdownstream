#!/usr/bin/env python3

import os

os.environ["NUMBA_CACHE_DIR"] = "./tmp/numba"
os.environ["MPLCONFIGDIR"] = "./tmp/matplotlib"

import scanpy as sc
import seaborn as sns
import matplotlib
import matplotlib.pyplot as plt
import numpy as np
import scipy.stats as stat
import platform
import base64
import json
import yaml
from matplotlib.pyplot import rc_context


def get_thresholds(adata, metric: str, nmads: [int, int]):
    """ Determines thresholds for outliers based on specified N median absolute deviations

    Parameters:
    adata: anndata object
    metric [str]: metric to be assessed, e.g. 'total_counts'
    namds [int, int]: number of median absolute deviations for [lower bound, upper bound]
    
    Returns:
    [int, int]: [lower bound, upper bound] threshold 

   """
    mat = adata.obs[metric]
    
    lower = np.median(mat) - nmads[0] * stat.median_abs_deviation(mat)
    upper = np.median(mat) + nmads[1] * stat.median_abs_deviation(mat) 
    
    return [lower, upper]



adata = sc.read_h5ad("${h5ad}")
prefix = "${prefix}"

## (currently) fixed parameters
selected_metrics = ['n_genes_by_counts', 'total_counts', 'pct_counts_mt']

namds_metric = {
    'n_genes_by_counts': [2.5, 5],
    'total_counts': [2.5, 5],
    'pct_counts_mt': [2.5,5]
}

dpi_plots = 120



# check for gene name, symbol etc. -- TODO: add as a data preparation step before calling QC
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

# get thresholds    
thresholds_metric = { metric:get_thresholds(adata, metric, namds_metric[metric]) for metric in selected_metrics}

# store threshold in anndata
adata.uns['n_mads_obs_thresholds'] = {
    'thresholds': thresholds_metric,
    'criteria_n_mads': namds_metric
}

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

# histograms
fig_width = 20
fig_height = 10

hist_color = sns.color_palette()[0]
line_color = sns.color_palette()[3]

with rc_context({'figure.figsize': (fig_width, fig_height)}):
    fig, (ax1, ax2, ax3) = plt.subplots(3)
    plt.subplots_adjust(hspace=0.5)
    sns.histplot(adata.obs['n_genes_by_counts'], kde=True, color=hist_color, edgecolor='black', ax=ax1)
    axis_ymax = ax1.get_ylim()[1]
    ax1.vlines(thresholds_metric['n_genes_by_counts'][0], 0, axis_ymax, colors=line_color, linestyles='solid')
    
    ax1.vlines(thresholds_metric['n_genes_by_counts'][1], 0, axis_ymax, colors=line_color, linestyles='solid')
    
    
    sns.histplot(adata.obs['total_counts'], kde=True, color=hist_color, edgecolor='black', ax=ax2)
    axis_ymax = ax2.get_ylim()[1]
    ax2.vlines(thresholds_metric['total_counts'][0], 0, axis_ymax, colors=line_color, linestyles='solid')
    ax2.vlines(thresholds_metric['total_counts'][1], 0, axis_ymax, colors=line_color, linestyles='solid')
    
    sns.histplot(adata.obs['pct_counts_mt'], kde=True, color=hist_color, edgecolor='black', ax=ax3)
    axis_ymax = ax3.get_ylim()[1]
    ax3.vlines(thresholds_metric['pct_counts_mt'][1], 0, axis_ymax, colors=line_color, linestyles='solid')
    
    path_plt_hist = "${prefix}_qc_histograms.png"
    fig.savefig(path_plt_hist)

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

    with open(path_plt_hist, "rb") as f_plot:
        image_string = base64.b64encode(f_plot.read()).decode("utf-8")
        image_html += f'<figure><div class="mqc-custom-content-image"><img src="data:image/png;base64,{image_string}" /></div>'
        image_html +='<figcaption>The histograms show gene counts, total counts, and the percentage of mitochondrial gene expression per cell. '
        image_html += 'Red lines indicate automatically determined thresholds based on N median absolute deviations (MADs).'
        
        image_html += 'Automatically determined thresholds:<br>'
        for metric in ['n_genes_by_counts', 'total_counts', 'pct_counts_mt']: 
            lower = thresholds_metric[metric][0]
            upper = thresholds_metric[metric][1]
            image_html += f'{metric}: ({lower}, {upper})<br>'

        image_html += 'The thresholds were determined using the following number of MADs:<br>'
        for metric in ['n_genes_by_counts', 'total_counts', 'pct_counts_mt']: 
            lower = namds_metric[metric][0]
            upper = namds_metric[metric][1]
            image_html += f'{metric}: ({lower}, {upper})<br>'
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
        "scanpy": sc.__version__,
        "matplotlib": matplotlib.__version__,
    }
}



with open("versions.yml", "w") as f:
    yaml.dump(versions, f)
