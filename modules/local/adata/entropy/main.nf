process ADATA_ENTROPY {
    tag "$meta.id"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container 'docker.io/gcfntnu/scanpy:1.11.4'

    input:
    tuple val(meta), path(h5ad)
    val(group_col)
    val(entropy_col)

    output:
    tuple val(meta), path("${prefix}.h5ad"), emit: h5ad
    path "${prefix}.pkl"                   , emit: obs
    path "${prefix}.png"                   , emit: plots, optional: true
    path "${prefix}_mqc.json"              , emit: multiqc_files, optional: true
    path "versions.yml"                    , emit: versions

    script:
    prefix = task.ext.prefix ?: "${meta.id}_entropy"
    plot_basis = task.ext.plot_basis ?: null
    template 'entropy.py'

    stub:
    prefix = task.ext.prefix ?: "${meta.id}_entropy"
    plot_basis = task.ext.plot_basis ?: null
    """
    touch ${prefix}.h5ad
    touch ${prefix}.pkl

    if [ ${plot_basis ? 'true' : 'false'} ]; then
        touch ${prefix}.png
        touch ${prefix}_mqc.json
    fi

    touch ${prefix}_mqc.json
    touch versions.yml
    """
}
