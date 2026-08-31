# UK DRI scdownstream: Documentation

The documentation is split into the following pages:

- [Usage](usage.md)
  - How the pipeline works, how to run each of its three entry points, the samplesheet and contrasts
    file formats, and a description of every parameter group.
- [Output](output.md)
  - The results each stage produces and how to interpret them.
- [UK DRI usage](ukdri.md)
  - Running the pipeline on the UK DRI HPC: pinned checkout, Apptainer image cache, SLURM job
    pattern and stage chaining.
- [Acknowledgements](../ACKNOWLEDGEMENTS.md)
  - Provenance of this fork and attribution to the original nf-core/scdownstream authors.

Parameters are defined in [`nextflow_schema.json`](../nextflow_schema.json); `nextflow run … --help`
lists them with their defaults.

Operational guidance for UK DRI users lives on the
[UK DRI Informatics wiki](https://wiki.informatics.ukdri.ac.uk/en/Pipelines/nfcore_scdownstream).
