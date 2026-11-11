process ADATA_UPSETGENES {
    tag "${meta.id}"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    // Built from containers/scdownstream-tools/Dockerfile.
    container 'docker.io/nhecker/scdownstream-tools:1.11.4-1'

    input:
    tuple val(meta), val(names), path(h5ads)

    output:
    tuple val(meta), path("*.png"), emit: plot, optional: true
    path ("*_mqc.json"), emit: multiqc_files, optional: true
    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    prefix = task.ext.prefix ?: "${meta.id}"
    template('upsetplot.py')

    stub:
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.png
    touch ${prefix}_mqc.json
    touch versions.yml
    """
}
