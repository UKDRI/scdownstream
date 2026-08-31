# Running UK DRI scdownstream at UK DRI

This page summarises how the pipeline is run on the UK DRI HPC. It is a pointer document — the
operational source of truth is the
[UK DRI Informatics wiki page](https://wiki.informatics.ukdri.ac.uk/en/Pipelines/nfcore_scdownstream).

For pipeline parameters and general usage, see [Usage](usage.md).

## What lives where

|                       | Path                                                |
| --------------------- | --------------------------------------------------- |
| Pipeline checkout     | `/nfsdata/scripts/nf-core/dev/scdownstream/main.nf` |
| Nextflow binary       | `/nfsdata/bin/nextflow-25.04.7-dist`                |
| Apptainer image cache | `/nfsdata/apptainer`                                |
| HCOP ortholog tables  | `/nfsdata/genome/hcop/`                             |
| Results               | somewhere under `/data/$USER/`                      |

Runs use the pinned checkout rather than `nextflow pull`, so everyone on a project is running the
same code.

## Environment

```bash
export NXF_SINGULARITY_CACHEDIR=/nfsdata/apptainer
export NXF_APPTAINER_CACHEDIR=/nfsdata/apptainer
export NXF_OPTS='-Xms1g -Xmx4g'
```

`--singularity_cache_dir` defaults to `$NXF_SINGULARITY_CACHEDIR`, so exporting it is enough to make
the pipeline resolve the custom `.sif` images locally instead of trying to pull them — see
[Local container images](usage.md#local-container-images).

## Invocation

Always `-profile apptainer,gpu`, always an explicit `-entry`, and one `params_<entry>.yml` per stage:

```bash
exec=/nfsdata/bin/nextflow-25.04.7-dist
main=/nfsdata/scripts/nf-core/dev/scdownstream/main.nf

$exec run $main \
    -entry qc_clustering \
    -profile apptainer,gpu \
    -params-file params_qc_clustering.yml \
    -resume
```

The `gpu` profile matters: scVI integration, CellBender and scAR all run substantially faster on a
GPU, and the `htc` partition provides them.

## SLURM

Nextflow runs as a driver job which submits the per-process jobs itself, so the `#SBATCH` header
sizes only the driver — not the work.

```bash
#SBATCH --partition=htc
#SBATCH --cpus-per-task=2
#SBATCH --time=48:00:00
```

Give the driver a generous `--time`: it must outlive every task it submits.

### Per-process resources

Per-process resources come from [`conf/base.config`](../conf/base.config), selected by process label.
To raise them, pass a config with `-c`:

```groovy title="custom.config"
process {
  withLabel:process_high {
    memory = 225.GB
  }
}
```

**Datasets above roughly 250,000 cells need this.** The merge, integration and clustering steps hold
the whole object in memory, and the default `process_high` allocation is not enough.

`--memory_scale` is a blunter alternative that multiplies every memory request in `base.config` by a
factor.

## Stage chaining

The three stages are sequential, each consuming the previous stage's `<name>_finalized.h5ad`:

```text
params_qc_clustering.yml       → results/qc_clustering/<name>_finalized.h5ad
params_downstream.yml          → results/downstream/<name>_finalized.h5ad
params_differential_genes.yml  → results/differential_genes/
```

Submit stage 2 only once stage 1's `.h5ad` exists. On SLURM you can queue them up front with a
dependency, since the output path is known in advance:

```bash
jid1=$(sbatch --parsable run_qc_clustering.sh)
jid2=$(sbatch --parsable --dependency=afterok:$jid1 run_downstream.sh)
sbatch --dependency=afterok:$jid2 run_differential_genes.sh
```

## Reference data

`--ortholog_hcop_directory` already defaults to `/nfsdata/genome/hcop/`, so LIANA+ works out of the
box on the cluster. CellTypist models are downloaded at runtime; if the compute nodes cannot reach the
internet, pass a local `.pkl` path instead. The same applies to celldex references, which accept
paths to pre-downloaded tar archives.

## Job script generation

Job scripts and validated `params_<entry>.yml` files can be generated with the
`nf-core_scdownstream` Claude skill in the UK DRI informatics pipeline skills repository, which
validates every parameter against this pipeline's `nextflow_schema.json` before writing the files.

> [!NOTE]
> The skill currently covers only `qc_clustering` and `downstream`. It has no
> `params_differential_genes.yml` or `run_..._differential_genes.sh` template, so stage 3 job scripts
> must be written by hand until the skill is updated.

## Off-cluster runs

Running outside UK DRI requires overriding the site-specific defaults:

- `--singularity_cache_dir` must point at a directory holding `scanpy_1.11.4_coreinf_0.3.sif`,
  `scanpy_1.11.4_coreinf_0.1.sif`, `decoupler_latest.sif` and `pydeseq2_latest.sif`. The report,
  enrichment and marker export processes have **no public registry fallback** — see
  [Local container images](usage.md#local-container-images).
- `--ortholog_hcop_directory` must point at your own [HCOP](https://www.genenames.org/tools/hcop/)
  tables.
- Drop `gpu` from `-profile` if no GPU is available.
