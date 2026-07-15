process DIFFERENTIAL_GENES_PER_CLUSTER {
    tag "${meta.id}"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container
        ? 'oras://community.wave.seqera.io/library/pydeseq2_anndata:ed0f23c01ea418d9'
        : 'community.wave.seqera.io/library/pydeseq2_anndata:ed0f23c01ea418d9'}"

    input:
    tuple val(meta), path(h5ad)
    val(cluster_col)
    val(condition_col)

    output:
    tuple val(meta), path("*.csv"), emit: csv
    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    prefix = task.ext.prefix ?: "${meta.id}"
    template('deseq2_diffgenes_per_cluster.py')
}

process DIFFERENTIAL_GENES_PER_CONTRAST {
    tag "${meta.id}"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container
        ? 'oras://community.wave.seqera.io/library/pydeseq2_anndata:ed0f23c01ea418d9'
        : 'community.wave.seqera.io/library/pydeseq2_anndata:ed0f23c01ea418d9'}"

    input:
    tuple val(meta), path(h5ad)
    val(condition_col)
    val(contrast)

    output:
    tuple val(meta), path("*.csv"), emit: csv
    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    prefix = task.ext.prefix ?: "${meta.id}"
    template('deseq2_diffgenes_per_contrast.py')
}

process DIFFERENTIAL_GENES_PER_CLUSTER_CONTRAST {
    tag "${meta.id}"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container
        ? 'oras://community.wave.seqera.io/library/pydeseq2_anndata:ed0f23c01ea418d9'
        : 'community.wave.seqera.io/library/pydeseq2_anndata:ed0f23c01ea418d9'}"

    input:
    tuple val(meta), path(h5ad)
    val(cluster_col)
    val(condition_col)

    output:
    tuple val(meta), path("*.csv"), emit: csv
    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    prefix = task.ext.prefix ?: "${meta.id}"
    template('deseq2_diffgenes_per_cluster_contrast.py')
}
