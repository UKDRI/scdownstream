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
min_auc = float("${min_auc}")
n_top = int("${n_top}")
pct_nz_group = float("${pct_nz_group}")
min_logfc = float("${min_logfc}")
# Optional minimum COSG score; empty (or null) means the filter is not used
min_cosg_score = "${min_cosg_score}".strip()
min_cosg_score = None if min_cosg_score in ("", "null") else float(min_cosg_score)
# Separator between the parts of a set id (genes<sep>cluster<sep><cluster>)
sep = "${sep}"

if sep == "":
    raise ValueError("The set id separator (--markers_sep) must not be empty")

METHOD_SUMMARIES = {
    "auc": (
        "scdownstream - Marker genes per cluster: Wilcoxon rank-sum test and AUC (presto statistics) "
        "comparing the cells of the cluster against all other cells; genes are ranked by AUC."
    ),
    "cosg": (
        "scdownstream - Marker genes per cluster are ranked by the COSG cosine-similarity score "
        "(cluster against all other cells), with the Wilcoxon rank-sum test and AUC (presto statistics) "
        "of the same comparison."
    ),
}


def group_table(rgg: dict, group: str) -> pd.DataFrame:
    """One group's results, read from the uns fields directly and kept in the stored ranking order.

    sc.get.rank_genes_groups_df would drop the auc and cosg_scores fields.
    """
    names = np.asarray(rgg["names"][group]).astype(str)
    table = pd.DataFrame({
        "names": names,
        "scores": np.asarray(rgg["scores"][group], dtype=np.float64),
        "auc": np.asarray(rgg["auc"][group], dtype=np.float64),
        "logfoldchanges": np.asarray(rgg["logfoldchanges"][group], dtype=np.float64),
        "pvals": np.asarray(rgg["pvals"][group], dtype=np.float64),
        "pvals_adj": np.asarray(rgg["pvals_adj"][group], dtype=np.float64),
    })
    if "cosg_scores" in rgg:
        table["cosg_scores"] = np.asarray(rgg["cosg_scores"][group], dtype=np.float64)
    table["pct_nz_group"] = rgg["pts"][group].reindex(names).to_numpy(dtype=np.float64)
    table["pct_nz_reference"] = rgg["pts_rest"][group].reindex(names).to_numpy(dtype=np.float64)
    return table


# main
if uns_key not in adata.uns:
    raise ValueError(f"adata.uns has no '{uns_key}' (--markers_uns_key)")
rgg = adata.uns[uns_key]
missing_fields = [field for field in ["names", "scores", "auc", "logfoldchanges", "pvals", "pvals_adj", "pts", "pts_rest"]
                  if field not in rgg]
if "auc" in missing_fields:
    raise ValueError(
        f"adata.uns['{uns_key}'] has no 'auc' field: it was not made by the presto marker step "
        f"(e.g. an object from an older pipeline version, which used sc.tl.rank_genes_groups). "
        f"Re-run --stage downstream to recompute the marker genes."
    )
if missing_fields:
    raise ValueError(f"adata.uns['{uns_key}'] lacks the fields {missing_fields}")
if min_cosg_score is not None and "cosg_scores" not in rgg:
    raise ValueError(
        "--markers_min_cosg_score needs --markers_method cosg; the marker results "
        f"(adata.uns['{uns_key}']) have no COSG scores"
    )

print(
    f"Filters: AUC >= {min_auc}, pct in >= {pct_nz_group}, log2FC >= {min_logfc}, COSG score: "
    + ("not used" if min_cosg_score is None else f">= {min_cosg_score}")
    + f"; top {n_top} per cluster"
)

obs_key = rgg['params']['groupby']
method = rgg['params']['method']
corr_method = rgg['params']['corr_method']
ranked_by = rgg['params'].get('ranked_by', 'auc')

params = {
    'method': method,
    'ranked_by': ranked_by,
    'corr_method': corr_method,
    'n_top': n_top,
    'min_auc': min_auc,
    'pct_nz_group': pct_nz_group,
    'min_logfc':  min_logfc,
    'min_cosg_score': min_cosg_score,
    'set_id_separator': sep
    }

cluster_markers = {}

for cluster in rgg['names'].dtype.names:

    # The id parts must not contain the separator, otherwise ids cannot be split again.
    if sep in str(cluster):
        raise ValueError(
            f"The {obs_key} value '{cluster}' contains the set id separator '{sep}', so its "
            f"set id could not be split; choose a different --markers_sep"
        )

    name_short = f"genes{sep}cluster{sep}{cluster}"
    name_long = f"Differential genes - {obs_key} - cluster {cluster}"

    # stored order is the ranking order (AUC, or the COSG score)
    df_markers = group_table(rgg, cluster)
    keep = ((df_markers['pct_nz_group'] >= pct_nz_group) &
            (df_markers['logfoldchanges'] >= min_logfc) &
            (df_markers['auc'] >= min_auc))
    if min_cosg_score is not None:
        keep &= df_markers['cosg_scores'] >= min_cosg_score
    df_markers = df_markers[keep]
    df_markers = df_markers.head(n=n_top)

    if df_markers.shape[0] > 0:

        cluster_markers[name_short] = {
            'description': name_long,
            'list_name': df_markers['names'].tolist(),
            'list_p_value': df_markers['pvals'].tolist(),
            'list_p_value_adjusted': df_markers['pvals_adj'].tolist(),
            'list_logFC': df_markers['logfoldchanges'].tolist(),
            'list_score': df_markers['scores'].tolist(),
            'list_auc': df_markers['auc'].tolist(),
            'list_pct_nz_group': df_markers['pct_nz_group'].tolist(),
            'list_pct_nz_reference': df_markers['pct_nz_reference'].tolist(),
        }
        if 'cosg_scores' in df_markers:
            cluster_markers[name_short]['list_cosg_score'] = df_markers['cosg_scores'].tolist()

marker_sets = {
    'project_name': project_name,
    'sets': cluster_markers,
    'params': params,
    'method_summary': METHOD_SUMMARIES.get(ranked_by, METHOD_SUMMARIES['auc'])
}

# write output json
with gzip.open(f"{prefix}.json.gz", 'wt') as file:
    json.dump(marker_sets, file)
