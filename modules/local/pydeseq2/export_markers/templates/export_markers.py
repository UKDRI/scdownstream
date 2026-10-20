#!/usr/bin/env python3

import csv
import gzip
import json
import math
import os
import platform

import numpy as np
import pandas as pd
import yaml

# parameters
manifest_path = "${manifest}"
prefix = "${prefix}"
project_name = "${project_name}"
group_col = "${group_col}"
thr_adj_pvalue = float("${thr_adj_pvalue}")
n_top = int("${n_top}")
min_logfc = float("${min_logfc}")

# Columns written by DIFFERENTIAL_GENES_PER_CONTRAST (pydeseq2 DeseqStats.results_df).
numeric_cols = ["baseMean", "log2FoldChange", "lfcSE", "stat", "pvalue", "padj"]
# Headerless manifest, one row per result file.
manifest_cols = ["file", "subset", "contrast", "variable", "target_group", "reference_group"]
# pandas to_csv writes NaN as an empty field; only those (and explicit NaN spellings) are
# missing values. Restricting NA parsing to the numeric columns keeps a gene called "NA"
# or "NULL" as a feature name instead of turning it into NaN.
na_strings = ["", "NaN", "nan", "NA"]


def to_float(x):
    """Python float, or None for NaN/inf so that the JSON stays valid."""
    if x is None:
        return None
    x = float(x)
    return x if math.isfinite(x) else None


def read_manifest(path):
    if os.path.getsize(path) == 0:
        return []
    rows = []
    with open(path, newline="") as f:
        for line_no, row in enumerate(csv.reader(f, delimiter="\\t"), start=1):
            if not row or all(field.strip() == "" for field in row):
                continue
            if len(row) != len(manifest_cols):
                raise ValueError(
                    f"Manifest line {line_no} has {len(row)} field(s), expected "
                    f"{len(manifest_cols)} ({', '.join(manifest_cols)}): {row}"
                )
            rows.append(dict(zip(manifest_cols, (field.strip() for field in row))))
    return rows


def read_results(path):
    df = pd.read_csv(
        path,
        sep="\\t",
        dtype={"feature": str},
        keep_default_na=False,
        na_values={col: na_strings for col in numeric_cols},
    )
    missing = [col for col in ["feature"] + numeric_cols if col not in df.columns]
    if missing:
        raise ValueError(f"{path} is missing column(s): {', '.join(missing)}")
    for col in numeric_cols:
        df[col] = pd.to_numeric(df[col], errors="coerce")
    # pydeseq2 independent filtering leaves padj as NaN for genes that were not tested.
    return df[df["padj"].notna()]


def to_set(df, description):
    return {
        "description": description,
        "list_name": [str(x) for x in df["feature"]],
        "list_p_value": [to_float(x) for x in df["pvalue"]],
        "list_p_value_adjusted": [to_float(x) for x in df["padj"]],
        "list_logFC": [to_float(x) for x in df["log2FoldChange"]],
        "list_score": [to_float(x) for x in df["stat"]],
        "list_base_mean": [to_float(x) for x in df["baseMean"]],
    }


# main
marker_sets_by_id = {}

for entry in read_manifest(manifest_path):
    subset = entry["subset"]
    contrast = entry["contrast"]
    variable = entry["variable"]
    target = entry["target_group"]
    reference = entry["reference_group"]

    df = read_results(entry["file"])
    significant = df[(df["padj"] < thr_adj_pvalue) & (df["log2FoldChange"].abs() >= min_logfc)]

    df_up = (
        significant[significant["log2FoldChange"] > 0]
        .sort_values("stat", ascending=False)
        .head(n_top)
    )
    df_down = (
        significant[significant["log2FoldChange"] < 0]
        .sort_values("stat", ascending=True)
        .head(n_top)
    )

    base_description = (
        f"Differential genes - {group_col} {subset} - {contrast}: "
        f"{variable} {target} vs {reference}"
    )

    for direction, df_set, label in [
        ("up", df_up, f"higher in {target}"),
        ("down", df_down, f"lower in {target}"),
    ]:
        # Empty sets are omitted, as in the downstream marker JSON.
        if df_set.shape[0] == 0:
            continue
        set_id = f"genes_{subset}_{contrast}_{direction}"
        if set_id in marker_sets_by_id:
            raise ValueError(f"Duplicate marker set id '{set_id}' (manifest file {entry['file']})")
        marker_sets_by_id[set_id] = to_set(df_set, f"{base_description} - {label}")

params = {
    "method": "pydeseq2 Wald test on pseudobulk profiles",
    "corr_method": "benjamini-hochberg",
    "n_top": n_top,
    "threshold_adj_pvalue": thr_adj_pvalue,
    "min_logfc": min_logfc,
    "group_col": group_col,
}

marker_sets = {
    "project_name": project_name,
    "sets": marker_sets_by_id,
    "params": params,
    "method_summary": (
        "scdownstream - Differential genes per contrast are detected with PyDESeq2 on "
        "pseudobulk profiles (per sample and group label), tested within each instance "
        "of the group label."
    ),
}

# write output json; allow_nan=False guarantees valid JSON (NaN/inf were mapped to None)
with gzip.open(f"{prefix}.json.gz", "wt") as f:
    json.dump(marker_sets, f, allow_nan=False)

versions = {
    "${task.process}": {
        "python": platform.python_version(),
        "pandas": pd.__version__,
        "numpy": np.__version__,
    }
}

with open("versions.yml", "w") as f:
    yaml.dump(versions, f)
