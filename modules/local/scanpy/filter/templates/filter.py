#!/usr/bin/env python3

import os
import platform

os.environ["MPLCONFIGDIR"] = "./tmp/mpl"
os.environ["NUMBA_CACHE_DIR"] = "./tmp/numba"

import yaml
import base64
import json
import scanpy as sc
import seaborn as sns
import matplotlib.pyplot as plt
import numpy as np
import scipy.stats as stat
from threadpoolctl import threadpool_limits
from matplotlib.pyplot import rc_context

def get_thresholds(adata, metric: str, nmads: list[int, int]):
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

    lower = 0 if upper < 0 else lower
    upper = 0 if upper < 0 else upper

    return [lower, upper]


threadpool_limits(int("${task.cpus}"))
sc.settings.n_jobs = int("${task.cpus}")


adata = sc.read_h5ad("${h5ad}")
prefix = "${prefix}"
symbol_col = "${symbol_col}"
mito_genes = "${mito_genes}"
automatic_filtering = "${automatic_cell_filtering}" == "true"

## (currently) fixed parameters
selected_metrics = ['n_genes_by_counts', 'total_counts', 'pct_counts_mt']

namds_metric = {
    'n_genes_by_counts': [2.5, 5],
    'total_counts': [2.5, 5],
    'pct_counts_mt': [2.5,5]
}

dpi_plots = 120

#  TODO: move custom mito genes to QC metric step
#
#symbols = adata.var_names if symbol_col == "index" else adata.var[symbol_col]
#
#
# if mito_genes:
#     with open(mito_genes, "r") as f:
#         mito_genes = {
#             line.strip().lower()
#             for line in f
#             if line.strip() and not line.startswith("#")
#         }
#     adata.var["mt"] = symbols.str.lower().isin(mito_genes)
# else:
#     adata.var["mt"] = symbols.str.lower().str.startswith("mt-")

# sc.pp.calculate_qc_metrics(
#     adata, qc_vars=["mt"], percent_top=None, log1p=False, inplace=True
# )


# thresholds cells
if automatic_filtering:

    # get automatic thresholds    
    thresholds_metric = { metric:get_thresholds(adata, metric, namds_metric[metric]) for metric in selected_metrics}

    # store threshold in anndata
    adata.uns['n_mads_obs_thresholds'] = {
            'thresholds': thresholds_metric,
            'criteria_n_mads': namds_metric
            }

else:
    max_pct_mt = int("${max_mito_percentage}")
    min_pct_mt = 0
    min_counts = int("${min_counts}")
    if "${max_counts}" != "false":
        max_counts = int("${max_counts}")
    else:   
        max_counts = np.max(adata.obs['total_counts']) + 1
    min_genes = int("${min_genes}")
    if "${max_genes}" != "false":
        max_genes = int("${max_genes}")
    else:   
        max_genes = np.max(adata.obs['n_genes_by_counts']) + 1

    adata = adata[adata.obs.pct_counts_mt < max_pct_mt, :].copy()
    sc.pp.filter_cells(adata, min_counts=min_counts)
    sc.pp.filter_cells(adata, min_genes=min_genes)

    thresholds_metric = {
        'n_genes_by_counts': [min_genes, max_genes],
        'total_counts': [min_counts, max_counts],
        'pct_counts_mt': [min_pct_mt, max_pct_mt]
    }


## plots

# threshold histograms (needs to be done before filtering)
sc.set_figure_params(dpi=dpi_plots, fontsize=18, dpi_save=dpi_plots)

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
    
    path_plt_hist = "${prefix}_hist_thresholds.png"
    fig.savefig(path_plt_hist)


## filtering

# filter cells
min_genes, max_genes = thresholds_metric['n_genes_by_counts']
min_counts, max_counts = thresholds_metric['total_counts']
_ , max_pct_mt = thresholds_metric['pct_counts_mt']

adata = adata[adata.obs.pct_counts_mt < max_pct_mt, :].copy()
sc.pp.filter_cells(adata, min_genes=min_genes)
sc.pp.filter_cells(adata, max_genes=max_genes)
sc.pp.filter_cells(adata, min_counts=min_counts)
sc.pp.filter_cells(adata, max_counts=max_counts)

# filter genes
sc.pp.filter_genes(adata, min_counts=int("${min_counts_gene}"))
sc.pp.filter_genes(adata, min_cells=int("${min_cells}"))

## output

# save anndata
adata.write_h5ad(f"{prefix}.h5ad")

# MultiQC
with open("${prefix}_mqc.json", "w") as f_json:
    
    image_html = ""

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
    "${task.process}": {"python": platform.python_version(), "scanpy": sc.__version__}
}

with open("versions.yml", "w") as f:
    yaml.dump(versions, f)
