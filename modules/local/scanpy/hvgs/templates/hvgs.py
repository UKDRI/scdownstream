#!/usr/bin/env python3

import importlib.util
import os
import platform
from threadpoolctl import threadpool_limits

os.environ["MPLCONFIGDIR"] = "./tmp/mpl"
os.environ["NUMBA_CACHE_DIR"] = "./tmp/numba"

import scanpy as sc
import yaml

threadpool_limits(int("${task.cpus}"))
sc.settings.n_jobs = int("${task.cpus}")

adata = sc.read_h5ad("${h5ad}")
prefix = "${prefix}"
n_hvgs = int("${n_hvgs}")
batch_key = "${batch_key}"
subset_to_hvgs = "${subset_to_hvgs}" == "true"
flavor = "${flavor}"
# Fraction of the genes in each loess fit of flavor 'seurat_v3' (ignored by the other flavours)
span = float("${span}")

# Layer with the raw counts, which the count-based flavours read (X if the layer is missing)
counts_layer = "counts"


select_hvgs = adata.n_vars > n_hvgs and n_hvgs >= 0

kwargs = {}

if batch_key:
    kwargs["batch_key"] = batch_key

# If an actual limit is provided, use it
# Otherwise, scanpy will automatically determine the number of highly variable genes
if select_hvgs and n_hvgs > 0:
    kwargs["n_top_genes"] = n_hvgs

if flavor == "seurat":
    # Always run this: it also stores the statistics (uns["hvg"], means, dispersions) that
    # sc.pl.highly_variable_genes needs in the reports
    sc.pp.highly_variable_genes(adata, **kwargs)

elif flavor in ("pearson_residuals", "seurat_v3"):
    if counts_layer in adata.layers:
        kwargs["layer"] = counts_layer
    else:
        print(f"SCANPY_HVGS: no '{counts_layer}' layer, selecting HVGs on X, which must hold raw counts")

    if flavor == "pearson_residuals":
        # Needs an explicit number of genes; with no selection every gene is kept below
        if select_hvgs:
            if "n_top_genes" not in kwargs:
                raise ValueError("SCANPY_HVGS: flavor 'pearson_residuals' needs a number of HVGs greater than 0.")
            sc.experimental.pp.highly_variable_genes(adata, flavor=flavor, **kwargs)
    else:
        # Without n_top_genes scanpy silently selects 2000 genes for seurat_v3, so '0 = automatic'
        # is refused. Without selection it still runs, for the statistics the reports plot.
        if select_hvgs and "n_top_genes" not in kwargs:
            raise ValueError("SCANPY_HVGS: flavor 'seurat_v3' needs a number of HVGs greater than 0.")
        if importlib.util.find_spec("skmisc") is None:
            raise ImportError(
                "SCANPY_HVGS: flavor 'seurat_v3' needs the scikit-misc package, which this container does not "
                "have. Use flavor 'pearson_residuals', or a container with scikit-misc."
            )

        # scikit-misc fits one loess per batch; it crashes the interpreter (segfault, no exception)
        # on a batch with a single cell, so such batches are refused before the call
        if batch_key:
            cells_per_batch = adata.obs[batch_key].value_counts(ascending=True)
            cells_per_batch = cells_per_batch[cells_per_batch > 0]
            tiny_batches = cells_per_batch[cells_per_batch < 2]
            if len(tiny_batches) > 0:
                tiny_list = ", ".join(f"{name} ({count} cell)" for name, count in tiny_batches.items())
                raise ValueError(
                    f"SCANPY_HVGS: flavor 'seurat_v3' needs at least 2 cells per batch ('{batch_key}'), "
                    f"but these batches have fewer: {tiny_list}. Remove or merge them in the samplesheet, "
                    "or use flavour 'seurat' (--hvg_flavor seurat, or --train_hvg_flavor seurat for --stage train_models)."
                )
            cells_text = "cells per batch (smallest first): " + ", ".join(
                f"{name}: {count}" for name, count in cells_per_batch.items()
            )
        else:
            cells_text = f"cells: {adata.n_obs}"

        try:
            sc.pp.highly_variable_genes(adata, flavor=flavor, span=span, **kwargs)
        except ValueError as err:
            raise ValueError(
                f"SCANPY_HVGS: the seurat_v3 loess fit failed ({err}), typically on small or sparse batches. "
                f"Span {span:g}; {cells_text}. "
                "Fix: raise the span (--hvg_span 1, or --train_hvg_span 1 for --stage train_models) "
                "or use flavour 'seurat' (--hvg_flavor seurat / --train_hvg_flavor seurat), and rerun with -resume."
            ) from err

        # Recorded for the reports; uns["hvg"] is set by scanpy, setdefault in case it is not
        adata.uns.setdefault("hvg", {})["span"] = span

else:
    raise ValueError(
        f"SCANPY_HVGS: unknown flavor '{flavor}'; expected 'seurat', 'pearson_residuals' or 'seurat_v3'."
    )

if not select_hvgs:
    # Not more genes than requested (or selection disabled): keep every gene
    adata.var["highly_variable"] = True

adata.var[["highly_variable"]].to_pickle(f"{prefix}.pkl")

if subset_to_hvgs:
    adata = adata[:, adata.var["highly_variable"]]

adata.write_h5ad(f"{prefix}.h5ad")

# Versions

versions = {
    "${task.process}": {
        "python": platform.python_version(),
        "scanpy": sc.__version__
    }
}

with open("versions.yml", "w") as f:
    yaml.dump(versions, f)
