process SCANPY_COMBAT {
    tag "${meta.id}"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container 'docker.io/gcfntnu/scanpy:1.11.4'

    input:
    tuple val(meta), path(h5ad)
    val batch_col

    output:
    tuple val(meta), path("*.h5ad"), emit: h5ad
    path "*.pkl", emit: obsm
    path "*.npy", emit: layers
    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    prefix = task.ext.prefix ?: "${meta.id}"
    template('combat.py')

    stub:
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.h5ad
    touch ${prefix}.pkl
    touch ${prefix}.npy
    touch versions.yml
    """
}
