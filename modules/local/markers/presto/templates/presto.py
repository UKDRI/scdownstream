#!/usr/bin/env python3

# Marker genes per group: Wilcoxon rank-sum test and AUC of each group against all other cells.
# Python re-implementation of the method of presto's wilcoxauc (Korsunsky et al. 2019,
# https://github.com/immunogenomics/presto); the results match R presto 1.1.0 (AUC, p-values,
# BH-adjusted p-values and expressing fractions verified to ~1e-14).
#
# Method:
#   - Ranks and ties: only the non-zero values of a gene are ranked (average ranks, values
#     compared in single precision); their ranks are shifted up by the number of zeros, the zeros
#     share one average rank, and the zero block always counts as one tie group.
#   - U statistic of a group with n1 cells, nnz1 of them non-zero, against the n2 other cells:
#       U = (n1 - nnz1) * (n_zero + 1) / 2 + R1(non-zeros) - n1 * (n1 + 1) / 2
#   - AUC = U / (n1 * n2)
#   - p-values from the normal approximation, with N = n1 + n2 and tie group sizes t:
#       z = U - n1 * n2 / 2;  z = z - sign(z) * 0.5                     (continuity correction)
#       sigma = sqrt(n1 * n2 * ((N^3 - N) - sum(t^3 - t)) / (12 * (N^2 - N)))   (tie-corrected)
#       z = z / sigma;  z = 0 where not finite;  p = 2 * pnorm(-abs(z))
#   - FDR: Benjamini-Hochberg separately for each group, across all genes.
#   - Cells without a group label are left out, and only groups with cells are tested.
# Only the sparse rank sums and the tie term are computed here (scipy's mannwhitneyu needs dense
# data and is much slower); the rest comes from libraries: the group means and non-zero counts
# from sc.get.aggregate, the normal CDF from scipy.special.ndtr, the ranks from
# scipy.stats.rankdata and Benjamini-Hochberg from scipy.stats.false_discovery_control.
# Differences from presto: the fold change is scanpy's (log2 of the ratio of expm1 means, see
# below) instead of presto's difference of log means, the fractions of expressing cells are
# fractions (scanpy's pts) instead of percentages, and negative values are ranked below the
# zeros (presto assumes non-negative values). The expm1 assumes the natural log of
# sc.pp.log1p, which is how the pipeline normalises (modules/local/scanpy/normalization).
#
# Results are stored in adata.uns["rank_genes_groups"] in the structure of sc.tl.rank_genes_groups,
# so sc.pl.rank_genes_groups_dotplot, sc.get.rank_genes_groups_df and
# sc.tl.filter_rank_genes_groups read them unchanged.

import os
import platform

os.environ["NUMBA_CACHE_DIR"] = "./tmp/numba"
os.environ["MPLCONFIGDIR"] = "./tmp/matplotlib"

import anndata as ad
import numpy as np
import pandas as pd
import scanpy as sc
import scipy
import yaml
from scipy import sparse
from scipy.special import ndtr
from scipy.stats import false_discovery_control, rankdata

from threadpoolctl import threadpool_limits
threadpool_limits(int("${task.cpus}"))

adata = sc.read_h5ad("${h5ad}")
prefix = "${prefix}"
obs_key = "${obs_key}"

# Results go where sc.tl.rank_genes_groups puts them, so the dot plot, the enrichment and the
# marker export (--markers_uns_key) find them
uns_key = "rank_genes_groups"
# Genes per block when X is dense: bounds the dense working copy to cells x block genes. Only
# affects memory use, not the results.
dense_block_genes = 256
# Pseudocount of scanpy's fold change, kept identical so --enrich_min_fold_change and the
# exported logFC mean the same as with sc.tl.rank_genes_groups
fold_change_pseudocount = 1e-9


# Groups
if obs_key not in adata.obs.columns:
    raise ValueError(
        f"The grouping '{obs_key}' is not a column of adata.obs; "
        f"available columns: {', '.join(map(str, adata.obs.columns))}"
    )

labels = adata.obs[obs_key]
labelled = labels.notna().to_numpy()
counts = labels.value_counts(dropna=True)
if isinstance(labels.dtype, pd.CategoricalDtype):
    # category order, without categories that have no cells (presto drops empty levels)
    groups = [group for group in labels.cat.categories if counts.get(group, 0) > 0]
else:
    try:
        groups = sorted(counts.index)
    except TypeError:
        groups = sorted(counts.index, key=str)
group_names = [str(group) for group in groups]

if len(set(group_names)) != len(group_names):
    raise ValueError(f"The values of adata.obs['{obs_key}'] are not unique when written as text: {group_names}")

print(f"Groups in adata.obs['{obs_key}'] (cells):")
for group, name in zip(groups, group_names):
    print(f"  {name}: {counts[group]}")
n_unlabelled = int((~labelled).sum())
if n_unlabelled > 0:
    print(f"{n_unlabelled} cells have no '{obs_key}' label and are left out of the tests.")

if len(groups) < 2:
    print(f"Skipping the marker genes: adata.obs['{obs_key}'] has fewer than 2 groups.")
else:
    n_groups = len(groups)
    n_genes = adata.n_vars
    codes = pd.Categorical(labels[labelled], categories=groups).codes.astype(np.intp)
    n_cells = codes.size
    group_size = np.bincount(codes, minlength=n_groups).astype(np.float64)

    X = adata.X
    if X is None:
        raise ValueError("adata.X is empty; the marker genes are computed on adata.X")
    if not labelled.all():
        X = X[labelled]

    # Per gene and group: mean and number of non-zero cells, over the labelled cells only (an
    # AnnData that shares X). count_nonzero tests X != 0, so stored zeros count as zeros.
    # aggregate drops empty categories and may order the groups differently, so its rows are
    # matched to the groups by name.
    adata_labelled = ad.AnnData(X=X, obs=adata.obs.loc[labelled, [obs_key]])
    aggregated = sc.get.aggregate(adata_labelled, by=obs_key, func=["mean", "count_nonzero"])
    del adata_labelled
    group_rows = aggregated.obs_names.get_indexer(group_names)
    if (group_rows < 0).any() or aggregated.n_obs != n_groups:
        raise ValueError(
            f"sc.get.aggregate returned the groups {list(aggregated.obs_names)}, expected {group_names}"
        )
    mean_in = np.asarray(aggregated.layers["mean"][group_rows], dtype=np.float64).T
    nonzero_in = np.asarray(aggregated.layers["count_nonzero"][group_rows], dtype=np.float64).T
    del aggregated

    # Per gene and group: rank sum; per gene: tie term
    rank_sum = np.zeros((n_genes, n_groups), dtype=np.float64)
    tie_sum = np.zeros(n_genes, dtype=np.float64)

    def accumulate(gene: int, values: np.ndarray, value_codes: np.ndarray) -> None:
        """Statistics of one gene from its non-zero values and their group codes."""
        if not np.isfinite(values).all():
            raise ValueError(f"adata.X has NaN or infinite values in gene '{adata.var_names[gene]}'")
        # Values are compared in single precision, as presto does, so values that are equal as
        # float32 are ties
        keys = values.astype(np.float32)
        n_zero = n_cells - values.size
        n_negative = np.count_nonzero(values < 0)
        ranks = rankdata(keys)  # average ranks 1..nnz among the non-zero values
        ranks[values > 0] += n_zero  # the zeros rank below the positive values
        zero_rank = n_negative + (n_zero + 1) / 2  # average rank of the zero block
        nnz_group = np.bincount(value_codes, minlength=n_groups)
        rank_sum[gene] = (np.bincount(value_codes, weights=ranks, minlength=n_groups)
                          + (group_size - nnz_group) * zero_rank)
        _, tie_sizes = np.unique(keys, return_counts=True)
        tie_sizes = tie_sizes.astype(np.float64)
        n_zero = float(n_zero)
        # tie groups of size 1 add 0, so counting every distinct value is the same as presto
        # counting the groups of size > 1; the zero block is always one group
        tie_sum[gene] = np.sum(tie_sizes ** 3 - tie_sizes) + (n_zero ** 3 - n_zero)

    if sparse.issparse(X):
        # One CSC copy (same size as X) gives the non-zeros of each gene; explicit zeros are dropped
        # so they are counted with the zeros
        X_csc = X.tocsc()
        if X_csc is X:
            X_csc = X_csc.copy()
        X_csc.sum_duplicates()
        X_csc.eliminate_zeros()
        for gene in range(n_genes):
            start, end = X_csc.indptr[gene], X_csc.indptr[gene + 1]
            accumulate(gene, X_csc.data[start:end], codes[X_csc.indices[start:end]])
        del X_csc
    else:
        X = np.asarray(X)
        for block_start in range(0, n_genes, dense_block_genes):
            block = np.asfortranarray(X[:, block_start:block_start + dense_block_genes])
            for offset in range(block.shape[1]):
                column = block[:, offset]
                nonzero = np.flatnonzero(column)
                accumulate(block_start + offset, column[nonzero], codes[nonzero])
    del X

    # Wilcoxon rank-sum test and AUC (formulas in the header)
    n = float(n_cells)
    n_rest = n - group_size
    n1n2 = group_size * n_rest
    ustat = rank_sum - group_size * (group_size + 1) / 2
    auc = ustat / n1n2

    z = ustat - 0.5 * n1n2
    z = z - np.sign(z) * 0.5
    tie_factor = ((n ** 3 - n) - tie_sum) / (12 * (n ** 2 - n))
    sigma = np.sqrt(n1n2[np.newaxis, :] * tie_factor[:, np.newaxis])
    with np.errstate(divide="ignore", invalid="ignore"):
        z = z / sigma
    z[~np.isfinite(z)] = 0
    pvals = 2 * ndtr(-np.abs(z))
    # Benjamini-Hochberg for each group (column) separately, as presto does
    pvals_adj = false_discovery_control(pvals, axis=0, method="bh")

    # Fold change as scanpy computes it: log2 of the ratio of the expm1 means in and out of the group.
    # The mean of all other cells is the size-weighted mean of the other groups.
    sum_in = mean_in * group_size
    mean_rest = (sum_in.sum(axis=1, keepdims=True) - sum_in) / n_rest
    logfoldchanges = np.log2(
        (np.expm1(mean_in) + fold_change_pseudocount) / (np.expm1(mean_rest) + fold_change_pseudocount)
    )

    # Fractions of cells with a non-zero value in the group and in all other cells (scanpy's pts)
    pts = nonzero_in / group_size
    pts_rest = (nonzero_in.sum(axis=1, keepdims=True) - nonzero_in) / n_rest

    # Each group sorted by AUC (descending), then adjusted p-value (ascending), then log2 fold change
    # (descending); remaining ties keep the gene order of adata.var_names (np.lexsort is stable, its
    # last key is the primary one). The fold change is compared as stored (float32), so the stored
    # values explain the order; MARKERS_COSG re-sorts with the same keys.
    order = np.column_stack([
        np.lexsort((-logfoldchanges[:, g].astype(np.float32), pvals_adj[:, g], -auc[:, g]))
        for g in range(n_groups)
    ])
    gene_names = adata.var_names.astype(str).to_numpy()
    names = np.broadcast_to(gene_names[:, np.newaxis], (n_genes, n_groups))

    def per_group_records(values: np.ndarray, dtype: str) -> np.recarray:
        """genes x groups matrix as scanpy's record array: one field per group, in ranking order."""
        frame = pd.DataFrame({name: values[order[:, g], g] for g, name in enumerate(group_names)})
        return frame.to_records(index=False, column_dtypes={name: dtype for name in group_names})

    adata.uns[uns_key] = {
        "params": {
            "groupby": obs_key,
            "reference": "rest",
            "method": "presto",
            "use_raw": False,
            # No "layer": None here (no None in uns): anndata >= 0.12 writes None as a 'null' element that steps with anndata < 0.12 (ADATA_TORDS, ADATA_SPLITCOL) cannot read
            "corr_method": "benjamini-hochberg",
            "ranked_by": "auc",
            # Marker methods that contributed; MARKERS_COSG sets "presto,cosg"
            "methods": "presto",
        },
        "names": per_group_records(names, "O"),
        "scores": per_group_records(auc, "float32"),
        "auc": per_group_records(auc, "float64"),
        "logfoldchanges": per_group_records(logfoldchanges, "float32"),
        "pvals": per_group_records(pvals, "float64"),
        "pvals_adj": per_group_records(pvals_adj, "float64"),
        "pts": pd.DataFrame(pts, index=adata.var_names.astype(str), columns=group_names),
        "pts_rest": pd.DataFrame(pts_rest, index=adata.var_names.astype(str), columns=group_names),
    }

    adata.write_h5ad(f"{prefix}.h5ad")

# Versions

versions = {
    "${task.process}": {
        "python": platform.python_version(),
        "scanpy": sc.__version__,
        "anndata": ad.__version__,
        "numpy": np.__version__,
        "scipy": scipy.__version__,
        "pandas": pd.__version__,
    }
}

with open("versions.yml", "w") as f:
    yaml.dump(versions, f)
