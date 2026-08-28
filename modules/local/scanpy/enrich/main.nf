process SCANPY_ENRICH {
    tag "${meta.id}"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${params.singularity_cache_dir
        ? params.singularity_cache_dir + '/scanpy_1.11.4_coreinf_0.1.sif'
        : '/data/nhecker/apptainer/images/scanpy_1.11.4_coreinf_0.1.sif'}"

    input:
    tuple val(meta), path(h5ad)
    val(uns_key)
    val(species)
    val(min_in_group_fraction)
    val(min_fold_change)
    val(max_out_group_fraction)


    output:
    tuple val(meta), path("*.h5ad"), emit: h5ad
    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    prefix = task.ext.prefix ?: "${meta.id}"
    template('enrich.py')
}
