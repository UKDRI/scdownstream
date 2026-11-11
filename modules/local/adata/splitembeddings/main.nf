process ADATA_SPLITEMBEDDINGS {
    tag "${meta.id}"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container 'docker.io/gcfntnu/scanpy:1.11.4'

    input:
    tuple val(meta), path(h5ad)
    val embeddings

    output:
    tuple val(meta), path("*.h5ad"), emit: h5ad
    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    prefix = task.ext.prefix ?: "${meta.id}"
    template('split_embeddings.py')
}
