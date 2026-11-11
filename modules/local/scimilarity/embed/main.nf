process SCIMILARITY_EMBED {
    tag "${meta.id}"
    label 'process_medium'
    label 'process_gpu'

    conda "${moduleDir}/environment.yml"
    // Built from containers/scdownstream-scimilarity/Dockerfile.
    container 'docker.io/nhecker/scdownstream-scimilarity:0.3.0-1'

    input:
    tuple val(meta), path(h5ad)
    tuple val(meta2), path(model)

    output:
    tuple val(meta), path("${prefix}.h5ad"), emit: h5ad
    path ("X_${prefix}.pkl"), emit: obsm
    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    prefix = task.ext.prefix ?: "${meta.id}"
    template('embed.py')
}
