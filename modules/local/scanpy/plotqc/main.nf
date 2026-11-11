process SCANPY_PLOTQC {
    tag "${meta.id}"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container 'docker.io/gcfntnu/scanpy:1.11.4'

    input:
    tuple val(meta), path(h5ad)

    output:
    tuple val(meta), path("*.png"), emit: plots
    path ("*_mqc.json"), emit: multiqc_files
    path "versions.yml", emit: versions
    tuple val(meta), path("${prefix}.h5ad"), emit: h5ad

    when:
    task.ext.when == null || task.ext.when

    script:
    prefix = task.ext.prefix ?: "${meta.id}"
    section_name = task.ext.section_name ?: "QC Plots"
    description = task.ext.description ?: "Quality control plots"
    template('plotqc.py')

    stub:
    prefix = task.ext.prefix ?: "${meta.id}"
    section_name = task.ext.section_name ?: "QC Plots"
    description = task.ext.description ?: "Quality control plots"
    """
    touch ${prefix}_total_counts_vs_n_genes_by_counts.png
    touch ${prefix}_qc_violin_plots.png
    touch ${prefix}_mqc.json
    touch versions.yml
    touch ${prefix}.h5ad
    """
}
