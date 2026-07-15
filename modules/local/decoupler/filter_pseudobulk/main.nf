process FILTER_PSEUDOBULK {
    tag "${meta.id}"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container
        ? 'oras://community.wave.seqera.io/library/decoupler_anndata:1f90c5f9ae5f4f8d'
        : 'community.wave.seqera.io/library/decoupler_anndata:1f90c5f9ae5f4f8d'}"

    input:
    tuple val(meta), path(h5ad)
    val(min_counts)
    val(min_cells)

    output:
    tuple val(meta), path("*.h5ad"), emit: h5ad
    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    prefix = task.ext.prefix ?: "${meta.id}"
    template('filter_pseudobulk.py')
}
