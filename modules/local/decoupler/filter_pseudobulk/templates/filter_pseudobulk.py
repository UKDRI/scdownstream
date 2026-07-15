#!/usr/bin/env python3

import platform
import warnings
import yaml

warnings.filterwarnings("ignore")

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
