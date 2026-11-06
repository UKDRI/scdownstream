#!/usr/bin/env python3

# Re-ranks the presto marker genes (MARKERS_PRESTO) by the COSG score (Dai et al., cosine
# similarity between a gene and a group indicator, penalised by the gene's similarity to the
# other groups; https://github.com/genecell/COSG). The presto statistics stay; the COSG score is
# added as "cosg_scores", becomes "scores", and every per-group field is re-sorted by it.
#
# The BioContainer has no scanpy and no PyYAML: read and write with anndata, and write
# versions.yml by hand.

import os
import platform
from importlib.metadata import version

os.environ["NUMBA_CACHE_DIR"] = "./tmp/numba"
os.environ["MPLCONFIGDIR"] = "./tmp/matplotlib"

import anndata as ad
import numpy as np
import pandas as pd
import cosg

from threadpoolctl import threadpool_limits
threadpool_limits(int("${task.cpus}"))


def format_yaml_like(data: dict, indent: int = 0) -> str:
    """Formats a dictionary to a YAML-like string.

    Args:
        data (dict): The dictionary to format.
        indent (int): The current indentation level.

    Returns:
        str: A string formatted as YAML.
    """
    yaml_str = ""
    for key, value in data.items():
        spaces = "  " * indent
        if isinstance(value, dict):
            yaml_str += f"{spaces}{key}:\\n{format_yaml_like(value, indent + 1)}"
        else:
            yaml_str += f"{spaces}{key}: {value}\\n"
    return yaml_str


adata = ad.read_h5ad("${h5ad}")
prefix = "${prefix}"
obs_key = "${obs_key}"
mu = float("${mu}")

# The presto results, written by MARKERS_PRESTO under scanpy's rank_genes_groups key
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
if uns_key not in adata.uns or "auc" not in adata.uns[uns_key]:
    raise ValueError(
        f"adata.uns['{uns_key}'] has no presto results (no 'auc' field); MARKERS_COSG needs the "
        f"output of MARKERS_PRESTO"
    )
rgg = adata.uns[uns_key]
if rgg["params"].get("groupby") != obs_key:
    raise ValueError(
        f"The presto results are grouped by '{rgg['params'].get('groupby')}', not by '{obs_key}'"
    )
if obs_key not in adata.obs.columns:
    raise ValueError(
        f"The grouping '{obs_key}' is not a column of adata.obs; "
        f"available columns: {', '.join(map(str, adata.obs.columns))}"
    )

var_names = adata.var_names.astype(str)
if not var_names.is_unique:
    raise ValueError("adata.var_names are not unique, so the COSG scores cannot be matched to the genes")

group_names = list(rgg["names"].dtype.names)

# COSG on the same cells and groups as presto: cells without a label are left out (COSG
# raises on them), and the labels are the presto group names, in presto's order
labels = adata.obs[obs_key]
labelled = labels.notna().to_numpy()
label_text = labels[labelled].astype(str).to_numpy()
if set(label_text) != set(group_names):
    raise ValueError(
        f"The groups of adata.obs['{obs_key}'] ({sorted(set(label_text))}) differ from the groups of "
        f"the presto results ({group_names})"
    )

# A light AnnData that only references adata.X (no copy when every cell has a label)
adata_cosg = ad.AnnData(
    X=adata.X if labelled.all() else adata.X[labelled],
    obs=pd.DataFrame(
        {obs_key: pd.Categorical(label_text, categories=group_names)},
        index=adata.obs_names[labelled],
    ),
    var=pd.DataFrame(index=var_names),
)

print(f"COSG on adata.X, grouped by '{obs_key}' ({len(group_names)} groups, {adata_cosg.n_obs} cells), mu = {mu}")

cosg.cosg(
    adata_cosg,
    groupby=obs_key,
    mu=mu,
    remove_lowly_expressed=remove_lowly_expressed,
    n_genes_user=adata_cosg.n_vars,
    key_added=cosg_key,
    use_raw=False,
    layer=None,
    calculate_logfoldchanges=False,
    calculate_pvalues=False,
    return_by_group=False,
)
cosg_result = adata_cosg.uns[cosg_key]

cosg_groups = [str(name) for name in cosg_result["names"].dtype.names]
if set(cosg_groups) != set(group_names):
    raise ValueError(f"COSG returned the groups {cosg_groups}, presto has {group_names}")

# Merge by gene name and re-sort every per-group field by the COSG score. Genes without a COSG
# score get 0, the score COSG itself gives a gene with no expression. Ties (e.g. all those zeros)
# keep the presto order, i.e. by AUC.
per_group_fields = [key for key, value in rgg.items() if is_per_group(value, group_names)]
reordered = {key: {} for key in per_group_fields}
cosg_sorted = {}
n_missing = 0

for cosg_name, name in zip(cosg_result["names"].dtype.names, cosg_groups):
    cosg_by_gene = pd.Series(
        np.asarray(cosg_result["scores"][cosg_name], dtype=np.float64),
        index=np.asarray(cosg_result["names"][cosg_name]).astype(str),
    )
    presto_names = np.asarray(rgg["names"][name]).astype(str)
    scores = cosg_by_gene.reindex(presto_names)
    n_missing += int(scores.isna().sum())
    scores = scores.fillna(0).to_numpy()

    order = np.argsort(-scores, kind="stable")
    for key in per_group_fields:
        reordered[key][name] = np.asarray(rgg[key][name])[order]
    cosg_sorted[name] = scores[order]

if n_missing > 0:
    print(f"{n_missing} gene-group pairs had no COSG score and got 0")

for key in per_group_fields:
    if key == "names":
        dtypes = {name: "O" for name in group_names}
        columns = {name: reordered[key][name].astype(str) for name in group_names}
    else:
        dtypes = {name: rgg[key].dtype[name] for name in group_names}
        columns = reordered[key]
    rgg[key] = pd.DataFrame(columns, columns=group_names).to_records(index=False, column_dtypes=dtypes)

cosg_records = pd.DataFrame(cosg_sorted, columns=group_names).to_records(
    index=False, column_dtypes={name: "float32" for name in group_names}
)
rgg["cosg_scores"] = cosg_records
rgg["scores"] = cosg_records.copy()
rgg["params"]["ranked_by"] = "cosg"
rgg["params"]["cosg_mu"] = mu

adata.uns[uns_key] = rgg
adata.uns.pop(cosg_key, None)

adata.write_h5ad(f"{prefix}.h5ad")

# Versions

versions = {
    "${task.process}": {
        "python": platform.python_version(),
        "cosg": version("cosg"),
        "anndata": ad.__version__,
        "numpy": np.__version__,
        "pandas": pd.__version__,
    }
}

with open("versions.yml", "w") as f:
    f.write(format_yaml_like(versions))
