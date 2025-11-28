#!/usr/bin/env python3

import os
import platform

os.environ["NUMBA_CACHE_DIR"] = "./tmp/numba"
os.environ["MPLCONFIGDIR"] = "./tmp/matplotlib"

import pandas as pd
import scanpy as sc
import liana as li

from threadpoolctl import threadpool_limits

threadpool_limits(int("${task.cpus}"))


def format_yaml_like(data: dict, indent: int = 0) -> str:
    """Formats a dictionary to a YAML-like string.

    Args:
        data (dict): The dictionary to format.
        indent (int): The current indentation level.

    Returns:
        str: A string formatted as YAML.
    """
    yaml_str = ""
    for key, value in data.items():
        spaces = "  " * indent
        if isinstance(value, dict):
            yaml_str += f"{spaces}{key}:\\n{format_yaml_like(value, indent + 1)}"
        else:
            yaml_str += f"{spaces}{key}: {value}\\n"
    return yaml_str


adata = sc.read_h5ad("${h5ad}")
prefix = "${prefix}"
obs_key = "${obs_key}"
species = "${species}"

dict_species = {
    'homo_sapiens': 'human',
    'homo sapiens': 'human',
    'hs': 'human',
    'mus_musculus': 'mouse',
    'mus musculus': 'mouse',
    'mm': 'mouse'
}

species = dict_species[species] if species in dict_species.keys() else species

# Getting ortholog mappings if needed
resource = None
if species != "human":
    try:
        resource = li.rs.select_resource('consensus')
        map_df = li.rs.get_hcop_orthologs(url=f'https://ftp.ebi.ac.uk/pub/databases/genenames/hcop/human_{species}_hcop_fifteen_column.txt.gz',
                                   columns=['human_symbol', 'mouse_symbol'],
                                   # NOTE: HCOP integrates multiple resource, so we can filter out mappings in at least 3 of them for confidence
                                   min_evidence=3
                                   )
        map_df = map_df.rename(columns={'human_symbol':'source', 'mouse_symbol':'target'})

        resource = li.rs.translate_resource(resource,
                                 map_df=map_df,
                                 columns=['ligand', 'receptor'],
                                 replace=True,
                                 # Here, we will be harsher and only keep mappings that don't map to more than 1 mouse gene
                                 one_to_many=1
                                 )
    except:
        print(f"WARNING. Failed to download or create HCOP orthlog mapping for {species}. Treating gene names as human gene names.")
        resource = None


if adata.obs[obs_key].nunique() > 1:
    #if (adata.X < 0).nnz == 0:
    #    sc.pp.log1p(adata)
    try:
        li.mt.rank_aggregate(
            adata, obs_key, use_raw=False, resource=resource, verbose=True, n_jobs=int("${task.cpus}")
        )
        df: pd.DataFrame = adata.uns["liana_res"]

        df.to_pickle(f"{prefix}.pkl")
        adata.write_h5ad(f"{prefix}.h5ad")

    except ValueError as e:
        if "cannot set a frame with no defined index and a scalar" in str(e):
            print(f"Error: {e}")
        else:
            raise e
else:
    print(
        f"Skipping rank aggregation because the column {obs_key} has only one unique value."
    )

# Versions

versions = {
    "python": platform.python_version(),
    "scanpy": sc.__version__,
    "liana": li.__version__,
}

with open("versions.yml", "w") as f:
    f.write(format_yaml_like(versions))
