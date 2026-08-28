#!/usr/bin/env python3

import platform
import re

import anndata as ad
import yaml

adata = ad.read_h5ad("${h5ad}")
column = "${column}"

assert column in adata.obs.columns, f"Column {column} not found in adata."


def sanitize(name):
    """Replace characters that are unsafe in a file name with underscores.

    Dots are replaced as well: downstream the subset name is recovered with
    `h5ad.simpleName`, which strips everything from the first dot onwards.
    """
    return re.sub(r"[^A-Za-z0-9_]", "_", str(name)) or "unnamed"


# Distinct obs values can sanitize to the same name (e.g. "T cell" and "T/cell"), which
# would silently overwrite a subset. Append a numeric suffix to keep every subset.
used = set()

for value in adata.obs[column].unique():
    adata_subset = adata[adata.obs[column] == value]
    name = sanitize(value)
    candidate = name
    suffix = 1
    while candidate in used:
        candidate = f"{name}_{suffix}"
        suffix += 1
    used.add(candidate)
    adata_subset.write_h5ad(f"{candidate}.h5ad")

# Versions

versions = {
    "${task.process}": {
        "python": platform.python_version(),
        "anndata": ad.__version__
    }
}

with open("versions.yml", "w") as f:
    yaml.dump(versions, f)
