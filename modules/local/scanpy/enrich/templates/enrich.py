#!/usr/bin/env python3

import os
import json
import platform
import base64

os.environ["NUMBA_CACHE_DIR"] = "./tmp/numba"
os.environ["MPLCONFIGDIR"] = "./tmp/matplotlib"

import numpy as np
import scanpy as sc
import pandas as pd
import yaml

from threadpoolctl import threadpool_limits
threadpool_limits(int("${task.cpus}"))
sc.settings.n_jobs = int("${task.cpus}")

# parameters
adata = sc.read_h5ad("${h5ad}")
prefix = "${prefix}"
species = "${species}"
uns_key = "${uns_key}"

min_in_group_fraction=float("${min_in_group_fraction}")
min_fold_change=float("${min_fold_change}")
max_out_group_fraction=float("${max_out_group_fraction}")

# TODO: add to parameters/command-line
uns_key_fil = uns_key + '_filtered'
# Adjusted p-value cutoff for the genes sent to g:Profiler; the default of sc.queries.enrich, which
# was called on the AnnData before, so the presto results are queried with the same genes
pval_cutoff = 0.05


dict_species = {
    'human': 'hsapiens',
    'homo_sapiens': 'hsapiens',
    'homo sapiens': 'hsapiens',
    'hs': 'hsapiens',
    'mouse': 'mmusculus',
    'mus_musculus': 'mmusculus',
    'mus musculus': 'mmusculus',
    'mm': 'mmusculus'
}
species = dict_species[species] if species in dict_species else species

# avoid 'nan' issue when calling filtering function 
adata_temp = adata.copy()

# filter genes
sc.tl.filter_rank_genes_groups(adata_temp, key=uns_key, key_added=uns_key_fil,
                              min_in_group_fraction=min_in_group_fraction, min_fold_change=min_fold_change, max_out_group_fraction=max_out_group_fraction)

# compute enrichment per group
dict_enrich = {}

# Gene list per group from the filtered results: the genes kept by the filter (non-NaN names), in the
# stored ranking order, and, when presto ran (pvals_adj exists), only those with pvals_adj < pval_cutoff.
# This is what sc.queries.enrich(adata, group, key=...) selected via sc.get.rank_genes_groups_df, which
# raises a KeyError on COSG-only results (no p-values).
filtered = adata_temp.uns[uns_key_fil]
has_pvals_adj = 'pvals_adj' in filtered

for grp in filtered['names'].dtype.names:

    names = pd.Series(filtered['names'][grp], dtype=object)
    keep = names.notna().to_numpy()
    if has_pvals_adj:
        keep &= np.asarray(filtered['pvals_adj'][grp], dtype=float) < pval_cutoff
    gene_list = names[keep].astype(str).tolist()

    try:
        enrichment = sc.queries.enrich(gene_list, org=species)
        enrichment['significant'] = [str(val) for val in enrichment['significant']]
        enrichment['parents'] = [",".join(golist) for golist in enrichment['parents']]
        dict_enrich[grp] = enrichment.copy()
    except:
        print(f"WARNING. 'enrich' query failed for '{grp}'.")
        dict_enrich[grp]  = pd.DataFrame({ 'query': ['empty set'] })

# store in uns
adata.uns[uns_key]['enrich'] = dict_enrich.copy()

# save anndata
adata.write_h5ad(f"{prefix}.h5ad")


# Versions

versions = {
    "${task.process}": {
        "python": platform.python_version(),
        "scanpy": sc.__version__,
    }
}

with open("versions.yml", "w") as f:
    yaml.dump(versions, f)
