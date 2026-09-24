process DIFFERENTIAL_GENES_PER_CLUSTER {
    tag "${meta.id}"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    // Built from the Dockerfile in this module directory.
    container 'nhecker/pydeseq2:0.1'

    input:
    tuple val(meta), path(h5ad)
    val(cluster_col)
    val(blocking_col)

    output:
    tuple val(meta), path("*.tsv"), emit: tsv
    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    prefix = task.ext.prefix ?: "${meta.id}"
    blocking_col = blocking_col ?: ''
    template('deseq2_diffgenes_per_cluster.py')

    stub:
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch "${prefix}_0.tsv"
    touch "${prefix}_1.tsv"
    touch "versions.yml"
    """
}

process DIFFERENTIAL_GENES_PER_CONTRAST {
    tag "${meta.id}"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    // Built from the Dockerfile in this module directory.
    container 'nhecker/pydeseq2:0.1'

    input:
    tuple val(meta), path(h5ad), val(contrast)
    val(min_samples)

    output:
    tuple val(meta), path("*.tsv"), emit: tsv, optional: true
    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    prefix = task.ext.prefix ?: "${meta.id}"
    contrast_name = contrast.contrast_name
    variable = contrast.variable
    target_group = contrast.target_group
    reference_group = contrast.reference_group
    blocked_variables = contrast.list_of_blocked_variables ?: ''
    template('deseq2_diffgenes_per_contrast.py')

    stub:
    prefix = task.ext.prefix ?: "${meta.id}"
    safe_name = contrast.contrast_name.toString().replaceAll('[~:+\\- /]', '_')
    """
    touch "${prefix}_${safe_name}.tsv"
    touch "versions.yml"
    """
}
