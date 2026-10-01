#!/usr/bin/env python3

import platform

import anndata as ad
import pandas as pd
import yaml

# The input is a symlink into the previous task's work directory; never write through it
if "${h5ad}" == "${prefix}.h5ad":
    raise ValueError("ADATA_ADDMETADATA: input and output are both named ${prefix}.h5ad; set a different ext.prefix")

metadata_file = "${metadata}"
sample_col = "${sample_col}"

# Check the raw header first: pandas would silently rename duplicate columns to "x.1".
# Read it as a data row (header=None) so the names are neither renamed nor turned into NaN.
header = pd.read_csv(
    metadata_file,
    sep="\\t",
    header=None,
    nrows=1,
    dtype=str,
    keep_default_na=False,
    encoding="utf-8-sig",
)
columns = [c.strip() for c in header.iloc[0].tolist()]

if "" in columns:
    raise ValueError(f"ADATA_ADDMETADATA: {metadata_file} has an empty column name in its header (a trailing tab?)")

duplicated_cols = sorted({c for c in columns if columns.count(c) > 1})
if duplicated_cols:
    raise ValueError(f"ADATA_ADDMETADATA: {metadata_file} has duplicate column names: {', '.join(duplicated_cols)}")

if sample_col not in columns:
    raise ValueError(
        f"ADATA_ADDMETADATA: sample column '{sample_col}' not found in {metadata_file}. "
        f"Columns: {', '.join(columns)}"
    )

# Infer column types, but always read the sample ids as strings.
# index_col=False stops a row with extra fields from shifting its first value into the index.
metadata = pd.read_csv(
    metadata_file,
    sep="\\t",
    header=0,
    names=columns,
    index_col=False,
    dtype={sample_col: str},
    encoding="utf-8-sig",
)

adata = ad.read_h5ad("${h5ad}")

if "sample" not in adata.obs:
    raise ValueError("ADATA_ADDMETADATA: obs has no 'sample' column; ADATA_ADDSAMPLE must run first")

# Sample ids must be present and unique
ids = metadata[sample_col].str.strip()
metadata[sample_col] = ids

empty_ids = ids.isna() | (ids == "")
if empty_ids.any():
    rows = (metadata.index[empty_ids.to_numpy()] + 1).tolist()
    raise ValueError(f"ADATA_ADDMETADATA: column '{sample_col}' has empty or NA values in data row(s) {rows}")

duplicated_ids = ids[ids.duplicated()].unique().tolist()
if duplicated_ids:
    raise ValueError(f"ADATA_ADDMETADATA: column '{sample_col}' has duplicate ids: {', '.join(duplicated_ids)}")

# Never overwrite existing obs columns
new_cols = [c for c in columns if c != sample_col]
if not new_cols:
    raise ValueError(f"ADATA_ADDMETADATA: {metadata_file} must have at least one column besides '{sample_col}'")

# Columns the pipeline creates itself later on (ADATA_UNIFY): "label" and "batch" would make it
# fail or silently change the batch key, "sample_original" is where it keeps an old sample column.
# "outlier" is written by SCANPY_FILTER with --filtering_keep_outliers.
reserved = [c for c in new_cols if c in ("batch", "label", "sample_original", "outlier")]
if reserved:
    raise ValueError(
        f"ADATA_ADDMETADATA: metadata column name(s) {', '.join(reserved)} are reserved by the pipeline "
        "(batch, label, sample_original, outlier); rename them in the TSV (e.g. 'batch' -> 'seq_batch')"
    )

clashing = [c for c in new_cols if c in adata.obs.columns]
if clashing:
    raise ValueError(
        f"ADATA_ADDMETADATA: these metadata columns already exist in obs: {', '.join(clashing)}; "
        "rename them in the TSV"
    )

# Every sample in the object needs a metadata row; unused rows are normal (one task per sample)
obs_sample = adata.obs["sample"].astype(str)
obs_ids = obs_sample.unique().tolist()
tsv_ids = set(ids)

missing = [s for s in obs_ids if s not in tsv_ids]
if missing:
    raise ValueError(
        f"ADATA_ADDMETADATA: sample(s) {', '.join(missing)} have no row in column '{sample_col}' of {metadata_file}"
    )

print(f"{len(tsv_ids) - len(obs_ids)} of {len(tsv_ids)} metadata rows not used (samples not in this object)")

# Map the per-sample rows onto cells
md = metadata.set_index(sample_col)
add = md.loc[obs_sample.values]
add.index = adata.obs_names

for col in new_cols:
    values = add[col]
    # Bool counts as numeric in pandas but not in ADATA_MERGE, so store it as a category too
    if pd.api.types.is_bool_dtype(values) or not pd.api.types.is_numeric_dtype(values):
        # String categories; missing values stay missing
        values = values.where(values.isna(), values.astype(str)).astype("category")
    adata.obs[col] = values.values

print(f"Sample(s): {', '.join(obs_ids)}")
print(f"Columns added: {', '.join(new_cols)}")

adata.write_h5ad("${prefix}.h5ad")

# Versions

versions = {
    "${task.process}": {
        "python": platform.python_version(),
        "anndata": ad.__version__,
        "pandas": pd.__version__,
    }
}

with open("versions.yml", "w") as f:
    yaml.dump(versions, f)
