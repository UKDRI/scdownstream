process DOUBLET_REMOVAL {
    tag "$meta.id"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    // Built from containers/scdownstream-tools/Dockerfile.
    container 'docker.io/nhecker/scdownstream-tools:1.11.4-1'

    input:
    tuple val(meta), path(h5ad), path(predictions)
    val(threshold)

    output:
    tuple val(meta), path("${prefix}.h5ad"), emit: h5ad
    path("*_mqc.json")             , emit: multiqc_files, optional: true
    path "versions.yml"            , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    prefix    = task.ext.prefix    ?: "${meta.id}"
    template 'doublet_removal.py'

    stub:
    prefix = task.ext.prefix    ?: "${meta.id}"
    """
    touch ${prefix}.h5ad
    touch ${prefix}_mqc.json
    touch versions.yml
    """
}
