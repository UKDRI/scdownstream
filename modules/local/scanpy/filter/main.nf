process SCANPY_FILTER {
    tag "${meta.id}"
    label 'process_low'

    conda "${moduleDir}/environment.yml"
    container 'docker.io/gcfntnu/scanpy:1.11.4'

    input:
    tuple val(meta), path(h5ad)
    val symbol_col
    val min_cells
    val min_genes
    val max_genes
    val min_counts
    val max_counts
    val min_counts_gene
    val max_mito_percentage
    val automatic_cell_filtering

    output:
    tuple val(meta), path("${prefix}.h5ad"), emit: h5ad
    path ("*_mqc.json"), emit: multiqc_files
    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    prefix = task.ext.prefix ?: "${meta.id}"
    if ("${prefix}.h5ad" == "${h5ad}") {
        error("Input and output names are the same, use \"task.ext.prefix\" to disambiguate!")
    }
    section_name = task.ext.section_name ?: "Filtering"
    description = task.ext.description ?: "Thresholds applied for filtering cells."
    template('filter.py')

    stub:
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.h5ad
    touch ${prefix}_filtering_thresholds.png
    touch ${prefix}_mqc.json
    touch versions.yml
    """
}
