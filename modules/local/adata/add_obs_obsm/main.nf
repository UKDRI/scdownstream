process ADATA_ADD_OBS_OBSM {
    tag "${meta.id}"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container
        ? 'oras://community.wave.seqera.io/library/anndata_pyyaml:5f82ece6392dc30c'
        : 'community.wave.seqera.io/library/anndata_pyyaml:b30e03a395613673'}"

    input:
    tuple val(meta), path(adata_base), path(adata_add)

    output:
    tuple val(meta), path("*.h5ad"), emit: h5ad
    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    prefix = task.ext.prefix ?: "${meta.id}"
    template('add_obs_obsm.py')
}