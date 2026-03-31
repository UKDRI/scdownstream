process ADATA_PUBLISH {
    tag "${meta.id}"
    label 'process_medium'

    input:
    tuple val(meta), path(h5ad)

    output:
    tuple val(meta), path("${prefix}.h5ad"), emit: h5ad

    when:
    task.ext.when == null || task.ext.when

    script:
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    mv ${h5ad} ${prefix}.h5ad
    """
}
