#!/usr/bin/env python3

import os
import platform
import warnings
import yaml

warnings.filterwarnings("ignore")

os.environ["NUMBA_CACHE_DIR"] = "./tmp/numba"
os.environ["MPLCONFIGDIR"] = "./tmp/matplotlib"
os.environ["XDG_CACHE_HOME"] = "./tmp/matplotlib/cache"

os.makedirs("./tmp/numba", exist_ok=True)
os.makedirs("./tmp/matplotlib/cache", exist_ok=True)

import anndata as ad
import decoupler as dc

adata = ad.read_h5ad("${h5ad}")
prefix = "${prefix}"
min_counts = int("${min_counts}")
min_cells = int("${min_cells}")

if not hasattr(adata, "obs"):
    raise ValueError("Input object is not an AnnData object")

dc.pp.filter_samples(
    adata,
    min_counts=min_counts,
    min_cells=min_cells,
)

adata.write_h5ad(f"{prefix}.h5ad")

versions = {
    "${task.process}": {
        "python": platform.python_version(),
        "anndata": ad.__version__,
        "decoupler": dc.__version__,
    }
}
with open("versions.yml", "w") as f:
    yaml.dump(versions, f)
