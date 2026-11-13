process SCDS {
    tag "${meta.id}"
    label 'process_medium'

    // Built from containers/seurat-ukdri-suite/Dockerfile.
    container 'docker.io/nhecker/seurat-ukdri-suite:5.5.1-1'

    input:
    tuple val(meta), path(h5ad)

    output:
    tuple val(meta), path("${prefix}.h5ad"), emit: h5ad
    tuple val(meta), path("${prefix}.csv"), emit: predictions
    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    prefix = task.ext.prefix ?: "${meta.id}"
    template('scds.R')

    stub:
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.h5ad
    touch ${prefix}.csv
    touch versions.yml
    """
}
