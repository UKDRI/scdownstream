process SCANPY_READH5 {
    tag "${meta.id}"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container 'docker.io/gcfntnu/scanpy:1.11.4'

    input:
    tuple val(meta), path(h5)

    output:
    tuple val(meta), path("*.h5ad"), emit: h5ad
    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    prefix = task.ext.prefix ?: "${meta.id}"
    template('readh5.py')

    stub:
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch "${prefix}.h5ad"
    touch "versions.yml"
    """
}
