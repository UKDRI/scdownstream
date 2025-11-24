#!/usr/bin/env python3

import os
import json
import platform
import base64

os.environ["NUMBA_CACHE_DIR"] = "./tmp/numba"
os.environ["MPLCONFIGDIR"] = "./tmp/matplotlib"

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

for grp in adata_temp.uns['rank_genes_groups']['names'].dtype.names:
    
    try:
        enrichment = sc.queries.enrich(adata_temp, grp, key=uns_key_fil, org=species)
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
