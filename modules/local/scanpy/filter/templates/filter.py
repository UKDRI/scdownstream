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

    # ignore NaN (pct_counts_mt of cells with zero counts); those cells still fail the filter mask
    median = np.nanmedian(mat)
    mad = stat.median_abs_deviation(mat, nan_policy="omit")

    lower = median - nmads[0] * mad
    upper = median + nmads[1] * mad

    lower = 0 if lower < 0 else lower
    upper = 0 if upper < 0 else upper

    return [lower, upper]


threadpool_limits(int("${task.cpus}"))
sc.settings.n_jobs = int("${task.cpus}")


adata = sc.read_h5ad("${h5ad}")
prefix = "${prefix}"
symbol_col = "${symbol_col}"
automatic_filtering = "${automatic_cell_filtering}" == "true"
keep_outliers = "${task.ext.keep_outliers ?: false}" == "true"

## (currently) fixed parameters
selected_metrics = ['n_genes_by_counts', 'total_counts', 'pct_counts_mt']

namds_metric = {
    'n_genes_by_counts': [2.5, 5],
    'total_counts': [2.5, 5],
    'pct_counts_mt': [2.5,5]
}

dpi_plots = 120


# thresholds cells
# metric -> (N-MAD lower bound, param name, floor) for lower bounds raised to a user floor
raised_floors = {}

if automatic_filtering:

    # get automatic thresholds
    n_mads_thresholds = {
        metric: [float(bound) for bound in get_thresholds(adata, metric, namds_metric[metric])]
        for metric in selected_metrics
    }

    # --min_genes / --min_counts act as floors for the N-MAD lower bounds: median - N*MAD
    # can be <= 0 (clamped to 0), which would let every barcode pass, including near-empty ones
    floors = {
        'n_genes_by_counts': int("${min_genes}"),
        'total_counts': int("${min_counts}")
    }
    floor_params = {'n_genes_by_counts': '--min_genes', 'total_counts': '--min_counts'}

    thresholds_metric = {metric: list(bounds) for metric, bounds in n_mads_thresholds.items()}
    for metric, floor in floors.items():
        n_mads_lower = n_mads_thresholds[metric][0]
        if floor > n_mads_lower:
            thresholds_metric[metric][0] = float(floor)
            raised_floors[metric] = (n_mads_lower, floor_params[metric], floor)
            print(f"{metric}: N-MAD lower bound {n_mads_lower} raised to {floor_params[metric]} {floor}")

    # store thresholds in anndata: 'thresholds' are the applied ones, 'n_mads_thresholds' before flooring
    adata.uns['n_mads_obs_thresholds'] = {
            'thresholds': thresholds_metric,
            'criteria_n_mads': namds_metric,
            'n_mads_thresholds': n_mads_thresholds,
            'floors': floors
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

# filter cells: one mask over all criteria; genes and counts per cell are counted on adata.X
min_genes, max_genes = thresholds_metric['n_genes_by_counts']
min_counts, max_counts = thresholds_metric['total_counts']
_ , max_pct_mt = thresholds_metric['pct_counts_mt']

pass_min_genes, n_genes = sc.pp.filter_cells(adata, min_genes=min_genes, inplace=False)
pass_max_genes, _ = sc.pp.filter_cells(adata, max_genes=max_genes, inplace=False)
pass_min_counts, n_counts = sc.pp.filter_cells(adata, min_counts=min_counts, inplace=False)
pass_max_counts, _ = sc.pp.filter_cells(adata, max_counts=max_counts, inplace=False)

cell_mask = (
    (adata.obs['pct_counts_mt'].to_numpy() < max_pct_mt)
    & pass_min_genes & pass_max_genes & pass_min_counts & pass_max_counts
)
n_outliers = int((~cell_mask).sum())

# the columns the in-place sc.pp.filter_cells calls add
adata.obs['n_genes'] = n_genes
adata.obs['n_counts'] = n_counts

if keep_outliers:
    if 'outlier' in adata.obs:
        print("WARNING: obs column 'outlier' already exists and is overwritten")
    adata.obs['outlier'] = ~cell_mask
else:
    adata = adata[cell_mask].copy()

adata.uns['cell_filtering'] = {
    'thresholds': {metric: [float(lower), float(upper)] for metric, (lower, upper) in thresholds_metric.items()},
    'keep_outliers': keep_outliers,
    'n_outliers': n_outliers
}

# filter genes
sc.pp.filter_genes(adata, min_counts=int("${min_counts_gene}"))
sc.pp.filter_genes(adata, min_cells=int("${min_cells}"))

# remove cells left without counts by the gene filter, regardless of keep_outliers: empty cells
# break downstream steps. This can drop genes below min_cells again, which is accepted (no iteration)
counts_per_cell = np.asarray(adata.X.sum(axis=1)).ravel()
empty_cells = counts_per_cell == 0
n_empty = int(empty_cells.sum())
print(f"{n_empty} cells removed with no counts left after gene filtering")
if n_empty > 0:
    adata = adata[~empty_cells].copy()
adata.uns['cell_filtering']['n_empty_after_gene_filtering'] = n_empty

# n_genes / n_counts of the remaining cells, counted on the final matrix
adata.obs['n_genes'] = np.asarray((adata.X > 0).sum(axis=1)).ravel()
adata.obs['n_counts'] = np.asarray(adata.X.sum(axis=1)).ravel()

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
        if automatic_filtering:
            image_html += 'Red lines indicate automatically determined thresholds based on N median absolute deviations (MADs).'
            for metric, (n_mads_lower, floor_param, floor) in raised_floors.items():
                image_html += f' Lower bounds for {metric} were raised from the N-MAD value {n_mads_lower} to the {floor_param} floor {floor}.'
        else:
            image_html += 'Red lines indicate default or user specified thresholds.'
        if keep_outliers:
            image_html += f' {n_outliers} cells marked as outliers (kept).'
        else:
            image_html += f' {n_outliers} cells removed.'
        image_html += f' {n_empty} cells removed with no counts left after gene filtering. '

        image_html += "Thresholds:<br>"
        image_html += "<table style='border-collapse: separate; border-spacing: 0.5em;'>"
        image_html += "<tr style='background-color: #F3F4F6'><th>Metric</th><th>lower bound</th><th>upper bound</th></tr>"
        for metric in ['n_genes_by_counts', 'total_counts', 'pct_counts_mt']: 
            lower = thresholds_metric[metric][0]
            upper = thresholds_metric[metric][1]
            image_html += f"<tr style='background-color: #E0F2FE'><td>{metric}</td><td>{lower}</td><td>{upper}</td></tr>"
        image_html += '</table><br>'
        
        if automatic_filtering:
            image_html += 'The thresholds were determined using the following number of MADs:<br>'
            image_html += "<table style='border-collapse: separate; border-spacing: 0.5em;'>"
            image_html += "<tr style='background-color: #F3F4F6'><th>Metric</th><th>lower bound</th><th>upper bound</th></tr>"
            for metric in ['n_genes_by_counts', 'total_counts', 'pct_counts_mt']: 
                lower = namds_metric[metric][0]
                upper = namds_metric[metric][1]
                image_html += f"<tr style='background-color: #E0F2FE'><td>{metric}</td><td>{lower}</td><td>{upper}</td></tr>"
            image_html += '</table><br>'
        
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
