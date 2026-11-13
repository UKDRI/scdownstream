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

PRESTO_TEXT = (
    "Wilcoxon rank-sum test and AUC (presto statistics) comparing the cells of the cluster against all "
    "other cells"
)
COSG_TEXT = "the COSG cosine-similarity score (cluster against all other cells)"
# One summary per (methods, ranked_by); the tie-break order is the one MARKERS_PRESTO / MARKERS_COSG sort by
METHOD_SUMMARIES = {
    ("presto", "auc"): (
        f"scdownstream - Marker genes per cluster: {PRESTO_TEXT}; genes are ranked by AUC, ties by adjusted "
        "p-value and then log2 fold change."
    ),
    ("presto,cosg", "auc"): (
        f"scdownstream - Marker genes per cluster: {PRESTO_TEXT}; genes are ranked by AUC, ties by adjusted "
        "p-value, log2 fold change and then the COSG score. The COSG cosine-similarity score of the same "
        "comparison is included."
    ),
    ("presto,cosg", "cosg"): (
        f"scdownstream - Marker genes per cluster are ranked by {COSG_TEXT}, ties by log2 fold change and "
        f"then the fraction of expressing cells in the cluster; the {PRESTO_TEXT} are included."
    ),
    ("cosg", "cosg"): (
        f"scdownstream - Marker genes per cluster are ranked by {COSG_TEXT}, ties by log2 fold change and "
        "then the fraction of expressing cells in the cluster. No Wilcoxon test or AUC was computed."
    ),
}


def group_table(rgg: dict, group: str, has_presto: bool, has_cosg: bool) -> pd.DataFrame:
    """One group's results, read from the uns fields directly and kept in the stored ranking order.

    Only the fields that exist: auc, pvals and pvals_adj with presto, cosg_scores with COSG.
    sc.get.rank_genes_groups_df would drop the auc and cosg_scores fields.
    """
    names = np.asarray(rgg["names"][group]).astype(str)
    fields = ["scores", "logfoldchanges"]
    if has_presto:
        fields += ["auc", "pvals", "pvals_adj"]
    if has_cosg:
        fields.append("cosg_scores")
    table = pd.DataFrame({"names": names})
    for field in fields:
        table[field] = np.asarray(rgg[field][group], dtype=np.float64)
    table["pct_nz_group"] = rgg["pts"][group].reindex(names).to_numpy(dtype=np.float64)
    table["pct_nz_reference"] = rgg["pts_rest"][group].reindex(names).to_numpy(dtype=np.float64)
    return table


# main
if uns_key not in adata.uns:
    raise ValueError(f"adata.uns has no '{uns_key}' (--markers_uns_key)")
rgg = adata.uns[uns_key]
# presto results (MARKERS_PRESTO) have 'auc', COSG results (MARKERS_COSG) 'cosg_scores'
has_presto = "auc" in rgg
has_cosg = "cosg_scores" in rgg
if not has_presto and not has_cosg:
    raise ValueError(
        f"adata.uns['{uns_key}'] has neither presto ('auc') nor COSG ('cosg_scores') results: it was not made "
        f"by the marker steps of this pipeline (e.g. an object from an older pipeline version, which used "
        f"sc.tl.rank_genes_groups). Re-run --stage downstream to recompute the marker genes."
    )
required_fields = ["names", "scores", "logfoldchanges", "pts", "pts_rest"]
if has_presto:
    required_fields += ["pvals", "pvals_adj"]
missing_fields = [field for field in required_fields if field not in rgg]
if missing_fields:
    raise ValueError(f"adata.uns['{uns_key}'] lacks the fields {missing_fields}")
if min_cosg_score is not None and not has_cosg:
    raise ValueError(
        "--markers_min_cosg_score is set, but the marker results "
        f"(adata.uns['{uns_key}']) have no COSG scores (COSG not in --markers_method, or an object from an "
        "older pipeline version). Add cosg to --markers_method, or unset --markers_min_cosg_score."
    )

# --markers_min_auc needs the AUC, which only presto computes
if not has_presto:
    print(f"Note: --markers_min_auc ({min_auc}) is not applied: presto did not run, so there is no AUC.")
    min_auc = None

print(
    "Filters: AUC: " + ("not used" if min_auc is None else f">= {min_auc}")
    + f", pct in >= {pct_nz_group}, log2FC >= {min_logfc}, COSG score: "
    + ("not used" if min_cosg_score is None else f">= {min_cosg_score}")
    + f"; top {n_top} per cluster"
)

obs_key = rgg['params']['groupby']
method = rgg['params']['method']
# Only presto has p-values and thus a multiple-testing correction
corr_method = rgg['params'].get('corr_method') if has_presto else None
ranked_by = rgg['params'].get('ranked_by', 'auc' if has_presto else 'cosg')
# Objects from before 'methods' was stored: derive it from the fields
methods = rgg['params'].get('methods', ",".join(m for m, used in (("presto", has_presto), ("cosg", has_cosg)) if used))

params = {
    'method': method,
    'methods': methods,
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

    # stored order is the ranking order (see METHOD_SUMMARIES); filters only on metrics that exist
    df_markers = group_table(rgg, cluster, has_presto, has_cosg)
    keep = ((df_markers['pct_nz_group'] >= pct_nz_group) &
            (df_markers['logfoldchanges'] >= min_logfc))
    if min_auc is not None:
        keep &= df_markers['auc'] >= min_auc
    if min_cosg_score is not None:
        keep &= df_markers['cosg_scores'] >= min_cosg_score
    df_markers = df_markers[keep]
    df_markers = df_markers.head(n=n_top)

    if df_markers.shape[0] > 0:

        marker_set = {
            'description': name_long,
            'list_name': df_markers['names'].tolist(),
        }
        if has_presto:
            marker_set['list_p_value'] = df_markers['pvals'].tolist()
            marker_set['list_p_value_adjusted'] = df_markers['pvals_adj'].tolist()
        marker_set['list_logFC'] = df_markers['logfoldchanges'].tolist()
        # the primary metric: AUC or COSG score (ranked_by)
        marker_set['list_score'] = df_markers['scores'].tolist()
        if has_presto:
            marker_set['list_auc'] = df_markers['auc'].tolist()
        if has_cosg:
            marker_set['list_cosg_score'] = df_markers['cosg_scores'].tolist()
        marker_set['list_pct_nz_group'] = df_markers['pct_nz_group'].tolist()
        marker_set['list_pct_nz_reference'] = df_markers['pct_nz_reference'].tolist()
        cluster_markers[name_short] = marker_set

marker_sets = {
    'project_name': project_name,
    'sets': cluster_markers,
    'params': params,
    'method_summary': METHOD_SUMMARIES.get(
        (methods, ranked_by),
        f"scdownstream - Marker genes per cluster from {methods}, ranked by {ranked_by}."
    )
}

# write output json
with gzip.open(f"{prefix}.json.gz", 'wt') as file:
    json.dump(marker_sets, file)
