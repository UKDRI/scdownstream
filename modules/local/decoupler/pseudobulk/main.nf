process DECOUPLER_PSEUDOBULK {
    tag "${meta.id}"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    // Public base image, no derived image needed: it already ships this template's full
    // dependency set (decoupler 2.1.1 — required for the dc.pp.* API used below — plus
    // anndata 0.12.2, pandas and pyyaml). Under singularity/apptainer Nextflow converts
    // this to a .sif in NXF_SINGULARITY_CACHEDIR on first use.
    container 'gcfntnu/scanpy:1.11.4'

    input:
    tuple val(meta), path(h5ad)
    val(sample_col)
    val(group_cols)

    output:
    tuple val(meta), path("*.h5ad"), emit: h5ad
    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    prefix = task.ext.prefix ?: "${meta.id}"
    template('pseudobulk.py')

    stub:
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch "${prefix}.h5ad"
    touch "versions.yml"
    """
}
