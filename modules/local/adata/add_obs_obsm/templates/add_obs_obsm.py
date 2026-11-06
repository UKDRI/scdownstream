#!/usr/bin/env python3

import platform
import sys
from pathlib import Path

import anndata as ad
import pandas as pd
import numpy as np
import yaml

adata = ad.read_h5ad("${adata_base}")
prefix = "${prefix}"
# Comma-separated prefixes of the obs columns to copy; empty: all columns
obs_prefixes = "${obs_prefixes}"
# Comma-separated obsm keys to copy, used only when filter_obsm is set (empty then means none)
filter_obsm = "${filter_obsm}" == "true"
obsm_keys = "${obsm_keys}"
overwrite = "${overwrite}" == "true"
remove_unused_categories = "${remove_unused_categories}" == "true"

# Folder the process stages the files to add into (stageAs: 'add/*')
add_dir = Path("add")
# Index name of the TSV of the copied obs columns
tsv_index_label = "barcode"


def warn(message):
    print(f"WARNING: ADATA_ADD_OBS_OBSM: {message}", file=sys.stderr)


def split_list(value):
    return [item.strip() for item in value.split(",") if item.strip()]


obs_prefix_list = tuple(split_list(obs_prefixes))
obsm_key_list = split_list(obsm_keys) if filter_obsm else None

add_files = sorted(str(path) for path in add_dir.iterdir() if path.is_file()) if add_dir.is_dir() else []
if not add_files:
    raise ValueError("ADATA_ADD_OBS_OBSM: no file to add from.")

# obs columns and obsm entries to copy from each file, with its cell names
to_add = []
for add_file in add_files:
    adata_add = ad.read_h5ad(add_file)
    add_cols = [col for col in adata_add.obs.columns if not obs_prefix_list or str(col).startswith(obs_prefix_list)]
    add_obs = adata_add.obs[add_cols].copy()
    if obsm_key_list is None:
        add_obsm = dict(adata_add.obsm)
    else:
        add_obsm = {key: adata_add.obsm[key] for key in obsm_key_list if key in adata_add.obsm}
    add_obsnames = [str(bc) for bc in adata_add.obs_names]
    del adata_add
    to_add.append((add_file, add_obs, add_obsm, add_obsnames))

if obsm_key_list:
    not_found = [key for key in obsm_key_list if not any(key in add_obsm for _, _, add_obsm, _ in to_add)]
    if not_found:
        warn(f"obsm key(s) {not_found} not found in any file to add from.")

# subset to the cells every file to add from has
keep = np.ones(adata.n_obs, dtype=bool)
for add_file, _, _, add_obsnames in to_add:
    in_add = adata.obs_names.isin(add_obsnames)
    if not in_add.all():
        warn(f"{int((~in_add).sum())} of {adata.n_obs} cells are not in {add_file} and are removed.")
    keep &= in_add
adata = adata[keep].copy()

# add to base adata; a column or key two files share is taken from the first file
added_obs = {}
added_obsm = {}
replaced = []
for add_file, add_obs, add_obsm, add_obsnames in to_add:
    dict_idx = {bc: idx for idx, bc in enumerate(add_obsnames)}
    order_add = [dict_idx[bc] for bc in adata.obs_names]

    rows = add_obs.iloc[order_add, :]
    for col in rows.columns:
        if col in added_obs:
            warn(f"obs column '{col}' is in {added_obs[col]} and {add_file}; keeping the one from {added_obs[col]}.")
            continue
        if col in adata.obs:
            if not overwrite:
                continue
            replaced.append(f"obs['{col}']")
        values = rows[col].copy()
        values.index = adata.obs_names
        adata.obs[col] = values
        added_obs[col] = add_file

    for emb, value in add_obsm.items():
        if emb in added_obsm:
            warn(f"obsm key '{emb}' is in {added_obsm[emb]} and {add_file}; keeping the one from {added_obsm[emb]}.")
            continue
        if emb in adata.obsm:
            if not overwrite:
                continue
            replaced.append(f"obsm['{emb}']")
        if isinstance(value, pd.DataFrame):
            value = value.iloc[order_add].copy()
            value.index = adata.obs_names
        else:
            value = value[order_add, :].copy()
        adata.obsm[emb] = value
        added_obsm[emb] = add_file

if remove_unused_categories:
    # e.g. a subset of clusters: no empty groups in the categorical columns
    trimmed = []
    for col in adata.obs.columns:
        if isinstance(adata.obs[col].dtype, pd.CategoricalDtype):
            categories = adata.obs[col].cat.categories
            adata.obs[col] = adata.obs[col].cat.remove_unused_categories()
            kept = categories.isin(adata.obs[col].cat.categories)
            if not kept.all():
                trimmed.append(col)
                # keep the colours stored per category in step with the categories
                colors_key = f"{col}_colors"
                if colors_key in adata.uns and len(adata.uns[colors_key]) == len(categories):
                    adata.uns[colors_key] = np.asarray(adata.uns[colors_key])[kept]
    if trimmed:
        print(f"Removed unused categories from {len(trimmed)} obs columns: {', '.join(map(str, trimmed))}")

if replaced:
    warn(f"replaced {len(replaced)} existing entries of the base object: {', '.join(replaced)}")
print(f"Added {len(added_obs)} obs columns and {len(added_obsm)} obsm keys to {adata.n_obs} cells.")

adata.write_h5ad(f"{prefix}.h5ad")

if obs_prefix_list:
    if not added_obs:
        warn(f"no obs column starting with {list(obs_prefix_list)} was copied.")
    adata.obs[list(added_obs)].to_csv(f"{prefix}_added_obs.tsv", sep="\\t", index_label=tsv_index_label)

# Versions

versions = {
    "${task.process}": {
        "python": platform.python_version(),
        "anndata": ad.__version__,
        "pandas": pd.__version__,
        "numpy": np.__version__
    }
}

with open("versions.yml", "w") as f:
    yaml.dump(versions, f)
