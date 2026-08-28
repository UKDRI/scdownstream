process SCANPY_EXPORT_MARKERS {
    tag "${meta.id}"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${params.singularity_cache_dir
        ? params.singularity_cache_dir + '/scanpy_1.11.4_coreinf_0.1.sif'
        : '/data/nhecker/apptainer/images/scanpy_1.11.4_coreinf_0.1.sif'}"

    input:
    tuple val(meta), path(h5ad)
    val(project_name)
    val(uns_key)
    val(thr_adj_pvalue)
    val(n_top)
    val(pct_nz_group)
    val(min_logfc)
    
    output:
    tuple val(meta), path("*.json.gz"), emit: json
    //path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    prefix = task.ext.prefix ?: "${meta.id}"
    template('export_markers.py')
}