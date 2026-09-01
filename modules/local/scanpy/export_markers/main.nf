process SCANPY_EXPORT_MARKERS {
    tag "${meta.id}"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    // Public base image, no derived image needed: this template only uses scanpy, pandas,
    // numpy, yaml and threadpoolctl, all of which the base already ships. Under
    // singularity/apptainer Nextflow converts it to a .sif in NXF_SINGULARITY_CACHEDIR.
    container 'gcfntnu/scanpy:1.11.4'

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