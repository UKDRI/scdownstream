process SCANPY_LOG_NORMALIZE {
    tag "${meta.id}"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container 'docker.io/gcfntnu/scanpy:1.11.4'

    input:
    tuple val(meta), path(h5ad)

    output:
    tuple val(meta), path("${prefix}.h5ad"), emit: h5ad
    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    prefix = task.ext.prefix ?: "${meta.id}"
    batch_key = task.ext.batch_key ?: ""
    template('normalize.py')

    stub:
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.h5ad
    touch versions.yml
    """
}