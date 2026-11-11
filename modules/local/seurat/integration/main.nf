process SEURAT_INTEGRATION {
    tag "${meta.id}"
    label 'process_medium'

    // Built from containers/seurat-ukdri-suite/Dockerfile.
    container 'docker.io/nhecker/seurat-ukdri-suite:5.5.1-1'

    input:
    tuple val(meta), path(h5ad)
    val(batch_col)

    output:
    tuple val(meta), path("${prefix}.h5ad"), emit: h5ad
    path "versions.yml", emit: versions

    script:
    prefix = task.ext.prefix ?: "${meta.id}"
    template('integration.R')

    stub:
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch "${prefix}.h5ad"
    touch "versions.yml"
    """
}
