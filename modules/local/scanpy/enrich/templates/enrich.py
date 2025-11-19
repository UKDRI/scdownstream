#!/usr/bin/env python3

import os
import json
import platform
import base64

os.environ["NUMBA_CACHE_DIR"] = "./tmp/numba"
os.environ["MPLCONFIGDIR"] = "./tmp/matplotlib"

import scanpy as sc
import yaml

from threadpoolctl import threadpool_limits
threadpool_limits(int("${task.cpus}"))
sc.settings.n_jobs = int("${task.cpus}")

# parameters
adata = sc.read_h5ad("${h5ad}")
prefix = "${prefix}"
species = "${species}"
uns_key = "${uns_key}"

# TODO: add to parameters/command-line
uns_key_fil = uns_key + '_filtered'

min_in_group_fraction=0.25
min_fold_change=1
max_out_group_fraction=0.5

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
    dict_enrich[grp] = sc.queries.enrich(adata_temp, grp, key=uns_key_fil, org=species)
    dict_enrich[grp]['significant'] = [str(val) for val in dict_enrich[grp]['significant']]
    dict_enrich[grp]['parents'] = [",".join(golist) for golist in dict_enrich[grp]['parents']]

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
