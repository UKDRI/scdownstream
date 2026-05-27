process SCANPY_FILTER {
    tag "${meta.id}"
    label 'process_low'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container
        ? 'oras://community.wave.seqera.io/library/pyyaml_scanpy:158b12038812cf13'
        : 'community.wave.seqera.io/library/pyyaml_scanpy:61c9ab8e312bbe0a'}"

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
    path mito_genes

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
