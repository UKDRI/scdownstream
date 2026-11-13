process ADATA_ADD_OBS_OBSM {
    tag "${meta.id}"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container 'docker.io/gcfntnu/scanpy:1.11.4'

    input:
    // adata_add: one or more files whose obs columns and obsm keys are added to adata_base
    tuple val(meta), path(adata_base), path(adata_add, stageAs: 'add/*')

    output:
    tuple val(meta), path("${prefix}.h5ad"), emit: h5ad
    // The copied obs columns, written only when ext.obs_prefixes is set
    tuple val(meta), path("${prefix}_added_obs.tsv"), emit: tsv, optional: true
    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    prefix = task.ext.prefix ?: "${meta.id}"
    // Comma-separated prefixes of the obs columns to copy; null: all columns
    obs_prefixes = task.ext.obs_prefixes ?: ''
    // Comma-separated obsm keys to copy; null: all keys, an empty string: none
    filter_obsm = task.ext.obsm_keys != null
    obsm_keys = task.ext.obsm_keys ?: ''
    // Replace obs columns and obsm keys the base already has (default: keep the base's)
    overwrite = task.ext.overwrite ?: false
    // Drop the categories no remaining cell has from every categorical obs column (default: keep them)
    remove_unused_categories = task.ext.remove_unused_categories ?: false

    if ("${adata_base}" == "${prefix}.h5ad") {
        error("Input and output names are the same, set prefix in module configuration to disambiguate!")
    }
    template('add_obs_obsm.py')

    stub:
    prefix = task.ext.prefix ?: "${meta.id}"
    obs_prefixes = task.ext.obs_prefixes ?: ''
    """
    touch ${prefix}.h5ad
    touch versions.yml

    if [ -n "${obs_prefixes}" ]; then
        touch ${prefix}_added_obs.tsv
    fi
    """
}
