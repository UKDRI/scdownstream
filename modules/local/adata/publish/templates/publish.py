#!/usr/bin/env python3

import platform

import anndata as ad
import h5py
import yaml

# The input is a symlink into the previous task's work directory; never write through it
if "${h5ad}" == "${prefix}.h5ad":
    raise ValueError("ADATA_PUBLISH: input and output are both named ${prefix}.h5ad; set a different ext.prefix")

# Write the final object gzip-compressed to save storage; it reads back like any h5ad
adata = ad.read_h5ad("${h5ad}")
adata.write_h5ad("${prefix}.h5ad", compression="gzip")

# Versions

versions = {
    "${task.process}": {
        "python": platform.python_version(),
        "anndata": ad.__version__,
        "h5py": h5py.__version__,
    }
}

with open("versions.yml", "w") as f:
    yaml.dump(versions, f)
