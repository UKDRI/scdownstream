#!/usr/bin/env python3

import platform

import anndata as ad
import h5py
import numpy as np
import scipy
import scipy.sparse as sp
import yaml

# The input is a symlink into the previous task's work directory; never write through it
if "${h5ad}" == "${prefix}.h5ad":
    raise ValueError("ADATA_PUBLISH: input and output are both named ${prefix}.h5ad; set a different ext.prefix")


def to_csr(x, name):
    """Return x as a csr_matrix if it is a numeric dense or non-CSR sparse matrix, else unchanged."""
    if x is None or isinstance(x, sp.csr_matrix):
        return x
    if not (isinstance(x, np.ndarray) or sp.issparse(x)) or x.dtype.kind not in "biufc":
        return x
    print(f"ADATA_PUBLISH: converted {name} from {type(x).__name__} to csr_matrix")
    return sp.csr_matrix(x)


adata = ad.read_h5ad("${h5ad}")

# Published objects are always sparse CSR, even if an upstream step (or the input) made them dense
adata.X = to_csr(adata.X, "X")
for key in list(adata.layers.keys()):
    adata.layers[key] = to_csr(adata.layers[key], f"layers['{key}']")
if adata.raw is not None:
    raw_X = to_csr(adata.raw.X, "raw.X")
    # raw.X has no setter; rebuild raw with the same var, varm and obs index
    if raw_X is not adata.raw.X:
        adata.raw = ad.AnnData(X=raw_X, var=adata.raw.var, varm=dict(adata.raw.varm), obs=adata.obs[[]])

# Write the final object gzip-compressed to save storage; it reads back like any h5ad
adata.write_h5ad("${prefix}.h5ad", compression="gzip")

# Versions

versions = {
    "${task.process}": {
        "python": platform.python_version(),
        "anndata": ad.__version__,
        "h5py": h5py.__version__,
        "scipy": scipy.__version__,
    }
}

with open("versions.yml", "w") as f:
    yaml.dump(versions, f)
