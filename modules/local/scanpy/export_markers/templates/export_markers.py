#!/usr/bin/env python3

import os
import json
import platform
import base64
import gzip

os.environ["NUMBA_CACHE_DIR"] = "./tmp/numba"
os.environ["MPLCONFIGDIR"] = "./tmp/matplotlib"

import numpy as np
import scanpy as sc
import pandas as pd
import yaml

from threadpoolctl import threadpool_limits
threadpool_limits(int("${task.cpus}"))
sc.settings.n_jobs = int("${task.cpus}")

# parameters
adata = sc.read_h5ad("${h5ad}")
prefix = "${prefix}"
project_name = "${project_name}"
uns_key = "${uns_key}"
thr_adj_pvalue = float("${thr_adj_pvalue}")
n_top = int("${n_top}")
pct_nz_group = float("${pct_nz_group}")
min_logfc = float("${min_logfc}")


# main
obs_key = adata.uns['rank_genes_groups']['params']['groupby']
method = adata.uns['rank_genes_groups']['params']['method']
corr_method = adata.uns['rank_genes_groups']['params']['corr_method']

params = { 
    'method': method,
    'corr_method': corr_method,
    'n_top': n_top,
    'threshold_adj_pvalue': thr_adj_pvalue,
    'pct_nz_group': pct_nz_group,
    'min_logfc':  min_logfc
    }

cluster_markers = {}

for cluster in np.unique(adata.obs[obs_key]):

    name_short = 'markers_cluster_' + cluster
    name_long= 'Marker genes - ' + obs_key + ' - cluster ' + cluster

    df_markers = sc.get.rank_genes_groups_df(adata, cluster, key=uns_key)
    df_markers.sort_values('scores', ascending=False, inplace=True)
    df_markers = df_markers[ (df_markers['pvals_adj'] < thr_adj_pvalue) & (df_markers['pct_nz_group'] >= pct_nz_group) &
                          (df_markers['logfoldchanges'] >= min_logfc) ]
    df_markers = df_markers.head(n=n_top)

    if df_markers.shape[0] > 0:

        cluster_markers[name_short] = {
            'description': name_long,
            'list_name': list(df_markers['names']),
            'list_p_value': list(df_markers['pvals']),
            'list_p_value_adjusted': list(df_markers['pvals_adj']),
            'list_logFC': list(df_markers['logfoldchanges']),
            'list_score': list(df_markers['scores']),
            'list_pct_nz_group':  list(df_markers['pct_nz_group']),
            'list_pct_nz_reference':  list(df_markers['pct_nz_reference']),
        }

marker_sets = {
    'project_name': project_name,
    'sets': cluster_markers,
    'params': params
}

# write output json
with gzip.open(f"{prefix}.json.gz", 'wt') as file:
    json.dump(marker_sets, file)
