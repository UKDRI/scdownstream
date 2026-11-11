#!/usr/bin/env python3

# Marker genes by the COSG score (Dai et al. 2022, cosine similarity between a gene and a group
# indicator, penalised by the gene's similarity to the other groups; https://github.com/genecell/COSG).
# Two cases, decided by the input presto_ran (set by the workflow):
# - presto_ran is true (MARKERS_PRESTO ran right before; its output is this input): the COSG score
#   is added to the presto results as "cosg_scores" (matched by gene name) and every per-group field
#   is re-sorted with the keys of rank_by. "scores" is the primary metric (AUC or COSG score).
# - presto_ran is false (--markers_method cosg): the results are built from COSG alone in the
#   structure of sc.tl.rank_genes_groups: names, scores and cosg_scores (the COSG score),
#   logfoldchanges (COSG's, computed like scanpy's), pts and pts_rest. No AUC and no p-values.
#   rank_by must be "cosg". Results already in the input (e.g. from an earlier downstream run) are
#   replaced.
# Ranking keys, in order (the stored order is the ranking used by the marker export and the report):
#   auc:  AUC desc, adjusted p-value asc, log2 fold change desc, COSG score desc
#   cosg: COSG score desc, log2 fold change desc, fraction of expressing cells in the group desc
# Remaining ties keep the gene order of adata.var_names.

import os
import platform
from collections.abc import Mapping
from importlib.metadata import version

os.environ["NUMBA_CACHE_DIR"] = "./tmp/numba"
os.environ["MPLCONFIGDIR"] = "./tmp/matplotlib"

import anndata as ad
import numpy as np
import pandas as pd
import scanpy as sc
import yaml
import cosg

from threadpoolctl import threadpool_limits
threadpool_limits(int("${task.cpus}"))


adata = sc.read_h5ad("${h5ad}")
prefix = "${prefix}"
obs_key = "${obs_key}"
mu = float("${mu}")
# Ranking of the stored results: "auc" (needs presto results) or "cosg"
rank_by = "${rank_by}".strip().lower()
# True when MARKERS_PRESTO ran right before and its output is this input
presto_ran = "${presto_ran}" == "true"

if rank_by not in ("auc", "cosg"):
    raise ValueError(f"rank_by must be 'auc' or 'cosg' (--markers_rank_by), not '{rank_by}'")

# Results go where sc.tl.rank_genes_groups (and MARKERS_PRESTO) put them, so the dot plot, the
# enrichment and the marker export find them
uns_key = "rank_genes_groups"
# Where COSG writes its own result; it is merged into uns_key and then dropped
cosg_key = "cosg"
# No expression filter in COSG, so every gene gets a score (COSG gives filtered genes -1). The
# marker export and the report filter on the fraction of expressing cells themselves.
remove_lowly_expressed = False


def is_per_group(value, group_names: list) -> bool:
    """True for scanpy's per-group record arrays (one field per group)."""
    return (isinstance(value, np.ndarray) and value.dtype.names is not None
            and set(value.dtype.names) == set(group_names))


# Checks
if obs_key not in adata.obs.columns:
    raise ValueError(
        f"The grouping '{obs_key}' is not a column of adata.obs; "
        f"available columns: {', '.join(map(str, adata.obs.columns))}"
    )

var_names = adata.var_names.astype(str)
if not var_names.is_unique:
    raise ValueError("adata.var_names are not unique, so the COSG scores cannot be matched to the genes")
n_genes = len(var_names)

# Groups as MARKERS_PRESTO defines them: the categories with cells, in category order (else the
# sorted values); cells without a label are left out (COSG raises on them)
labels = adata.obs[obs_key]
labelled = labels.notna().to_numpy()
counts = labels.value_counts(dropna=True)
if isinstance(labels.dtype, pd.CategoricalDtype):
    groups = [group for group in labels.cat.categories if counts.get(group, 0) > 0]
else:
    try:
        groups = sorted(counts.index)
    except TypeError:
        groups = sorted(counts.index, key=str)
group_names = [str(group) for group in groups]
n_groups = len(group_names)

if len(set(group_names)) != n_groups:
    raise ValueError(f"The values of adata.obs['{obs_key}'] are not unique when written as text: {group_names}")
if n_groups < 2:
    raise ValueError(f"adata.obs['{obs_key}'] has fewer than 2 groups with cells, so there are no marker genes")
n_unlabelled = int((~labelled).sum())
if n_unlabelled > 0:
    print(f"{n_unlabelled} cells have no '{obs_key}' label and are left out.")

codes = pd.Categorical(labels[labelled], categories=groups).codes.astype(np.intp)

X = adata.X
if X is None:
    raise ValueError("adata.X is empty; the COSG score is computed on adata.X")
if not labelled.all():
    X = X[labelled]

rgg = adata.uns.get(uns_key)

if presto_ran:
    # MARKERS_PRESTO ran right before, so its results must be in the input; anything else is a wiring bug
    params = rgg.get("params", {}) if isinstance(rgg, Mapping) else {}
    if params.get("method") != "presto" or params.get("groupby") != obs_key:
        raise ValueError(
            f"presto_ran is true, but adata.uns['{uns_key}'] holds no presto results grouped by '{obs_key}' "
            f"(method '{params.get('method')}', groupby '{params.get('groupby')}'). The input of MARKERS_COSG "
            f"must be the output of MARKERS_PRESTO when presto_ran is true."
        )
    if set(rgg["names"].dtype.names) != set(group_names):
        raise ValueError(
            f"The groups of the presto results {list(rgg['names'].dtype.names)} differ from the groups of "
            f"adata.obs['{obs_key}'] {group_names}"
        )
    print("Adding the COSG score to the presto results of MARKERS_PRESTO")
else:
    if rank_by != "cosg":
        raise ValueError(
            f"rank_by is '{rank_by}', but MARKERS_PRESTO did not run: ranking by AUC needs MARKERS_PRESTO "
            f"(--markers_method presto,cosg); with COSG alone, rank by 'cosg'"
        )
    if rgg is not None:
        method = rgg.get("params", {}).get("method") if isinstance(rgg, Mapping) else None
        print(f"Replacing the existing adata.uns['{uns_key}'] (method '{method}', e.g. from an earlier downstream run)")
    print("The marker results are built from COSG alone (no AUC, no p-values)")

# COSG on the labelled cells, with the group names as labels, in group order. A light AnnData that
# only references X; no uns, so COSG's fold change uses the natural log, as MARKERS_PRESTO does.
adata_cosg = ad.AnnData(
    X=X,
    obs=pd.DataFrame(
        {obs_key: pd.Categorical.from_codes(codes, categories=group_names)},
        index=adata.obs_names[labelled],
    ),
    var=pd.DataFrame(index=var_names),
)

# Fractions of cells with a non-zero value, genes x groups: in the group (pts) and in all other
# cells (pts_rest), as scanpy's pts / pts_rest and MARKERS_PRESTO. count_nonzero tests X != 0, so
# stored zeros count as zeros. aggregate drops empty categories and may order the groups
# differently, so its rows are matched to the groups by name. Taken before COSG runs on adata_cosg.
aggregated = sc.get.aggregate(adata_cosg, by=obs_key, func=["count_nonzero"])
group_rows = aggregated.obs_names.get_indexer(group_names)
if (group_rows < 0).any() or aggregated.n_obs != n_groups:
    raise ValueError(f"sc.get.aggregate returned the groups {list(aggregated.obs_names)}, expected {group_names}")
nonzero_in = np.asarray(aggregated.layers["count_nonzero"][group_rows], dtype=np.float64).T
del aggregated
group_size = np.bincount(codes, minlength=n_groups).astype(np.float64)
pts = nonzero_in / group_size
pts_rest = (nonzero_in.sum(axis=1, keepdims=True) - nonzero_in) / (codes.size - group_size)

print(f"COSG on adata.X, grouped by '{obs_key}' ({n_groups} groups, {adata_cosg.n_obs} cells), mu = {mu}")

cosg.cosg(
    adata_cosg,
    groupby=obs_key,
    mu=mu,
    remove_lowly_expressed=remove_lowly_expressed,
    n_genes_user=n_genes,
    key_added=cosg_key,
    use_raw=False,
    layer=None,
    # presto's fold change (the same formula) is kept when presto ran
    calculate_logfoldchanges=not presto_ran,
    calculate_pvalues=False,
    return_by_group=False,
)
cosg_result = adata_cosg.uns[cosg_key]
del adata_cosg, X

cosg_fields = {str(field): field for field in cosg_result["names"].dtype.names}
if set(cosg_fields) != set(group_names):
    raise ValueError(f"COSG returned the groups {list(cosg_fields)}, expected {group_names}")
if not presto_ran and "logfoldchanges" not in cosg_result:
    raise ValueError("COSG returned no log fold changes (calculate_logfoldchanges=True)")


def cosg_in_var_order(key: str, group: str, fill: float) -> tuple[np.ndarray, int]:
    """COSG's per-group field `key` in adata.var_names order; genes COSG did not return get `fill`."""
    field = cosg_fields[group]
    values = pd.Series(
        np.asarray(cosg_result[key][field], dtype=np.float32),
        index=np.asarray(cosg_result["names"][field]).astype(str),
    ).reindex(var_names)
    n_missing = int(values.isna().sum())
    return values.fillna(fill).to_numpy(dtype=np.float32), n_missing


# Per-group fields as genes x groups arrays in adata.var_names order (dtype per field), so that
# sorting with a stable sort leaves the remaining ties in var_names order
fields = {}
dtypes = {}
n_missing = 0

if presto_ran:
    for key in [key for key, value in rgg.items() if is_per_group(value, group_names) and key != "names"]:
        dtypes[key] = rgg[key].dtype[group_names[0]]
        fields[key] = np.empty((n_genes, n_groups), dtype=dtypes[key])
    for g, name in enumerate(group_names):
        position = var_names.get_indexer(np.asarray(rgg["names"][name]).astype(str))
        if (position < 0).any() or np.unique(position).size != n_genes:
            raise ValueError(f"The presto results of group '{name}' do not hold each gene of adata.var_names once")
        for key in fields:
            fields[key][position, g] = np.asarray(rgg[key][name])
fields["cosg_scores"] = np.empty((n_genes, n_groups), dtype=np.float32)
dtypes["cosg_scores"] = np.float32
if not presto_ran:
    fields["logfoldchanges"] = np.empty((n_genes, n_groups), dtype=np.float32)
    dtypes["logfoldchanges"] = np.float32

for g, name in enumerate(group_names):
    # Genes without a COSG score get 0, the score COSG itself gives a gene with no expression
    fields["cosg_scores"][:, g], missing = cosg_in_var_order("scores", name, 0.0)
    n_missing += missing
    if not presto_ran:
        fields["logfoldchanges"][:, g], _ = cosg_in_var_order("logfoldchanges", name, np.nan)

if n_missing > 0:
    print(f"{n_missing} gene-group pairs had no COSG score and got 0")

# "scores" is the primary metric: the AUC (as MARKERS_PRESTO stores it, float32) or the COSG score
fields["scores"] = fields["auc"].astype(np.float32) if rank_by == "auc" else fields["cosg_scores"].copy()
dtypes["scores"] = np.float32

# Ranking per group (np.lexsort: stable, the last key is the primary one; descending = negated).
# The keys are compared as stored, so the stored values explain the order.
order = np.empty((n_genes, n_groups), dtype=np.intp)
for g in range(n_groups):
    if rank_by == "auc":
        keys = (-fields["cosg_scores"][:, g], -fields["logfoldchanges"][:, g].astype(np.float32),
                fields["pvals_adj"][:, g], -fields["auc"][:, g])
    else:
        keys = (-pts[:, g], -fields["logfoldchanges"][:, g].astype(np.float32), -fields["cosg_scores"][:, g])
    order[:, g] = np.lexsort(keys)


def per_group_records(values: np.ndarray, dtype) -> np.recarray:
    """genes x groups matrix as scanpy's record array: one field per group, in ranking order."""
    frame = pd.DataFrame({name: values[order[:, g], g] for g, name in enumerate(group_names)})
    return frame.to_records(index=False, column_dtypes={name: dtype for name in group_names})


gene_names = np.broadcast_to(var_names.to_numpy()[:, np.newaxis], (n_genes, n_groups))
records = {"names": per_group_records(gene_names, "O")}
records.update({key: per_group_records(fields[key], dtypes[key]) for key in fields})

if presto_ran:
    rgg.update(records)
    rgg["params"]["ranked_by"] = rank_by
    rgg["params"]["methods"] = "presto,cosg"
    rgg["params"]["cosg_mu"] = mu
else:
    rgg = {
        "params": {
            "groupby": obs_key,
            "reference": "rest",
            "use_raw": False,
            # No "layer": None here (no None in uns): anndata >= 0.12 writes None as a 'null' element that steps with anndata < 0.12 (ADATA_TORDS, ADATA_SPLITCOL) cannot read
            "method": "cosg",
            "ranked_by": "cosg",
            "methods": "cosg",
            "cosg_mu": mu,
        },
        **records,
        "pts": pd.DataFrame(pts, index=var_names, columns=group_names),
        "pts_rest": pd.DataFrame(pts_rest, index=var_names, columns=group_names),
    }

print(f"Ranked by {rank_by}; methods: {rgg['params']['methods']}")

adata.uns[uns_key] = rgg
adata.uns.pop(cosg_key, None)

adata.write_h5ad(f"{prefix}.h5ad")

# Versions

versions = {
    "${task.process}": {
        "python": platform.python_version(),
        "scanpy": sc.__version__,
        "anndata": ad.__version__,
        "numpy": np.__version__,
        "pandas": pd.__version__,
        "cosg": version("cosg"),
    }
}

with open("versions.yml", "w") as f:
    yaml.dump(versions, f)
