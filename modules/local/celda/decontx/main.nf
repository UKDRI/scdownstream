process CELDA_DECONTX {
    tag "$meta.id"
    label 'process_medium'

    // Built from containers/seurat-ukdri-suite/Dockerfile.
    container 'docker.io/nhecker/seurat-ukdri-suite:5.5.1-1'

    input:
    tuple val(meta), path(h5ad), path(raw)
    val(batch_col)
    val(input_layer)
    val(output_layer)

    output:
    tuple val(meta), path("*.h5ad"), emit: h5ad
    path "versions.yml"            , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    prefix = task.ext.prefix ?: "${meta.id}"
    template 'decontx.R'

    stub:
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.h5ad
    touch versions.yml
    """
}
