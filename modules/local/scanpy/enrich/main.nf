process SCANPY_ENRICH {
    tag "${meta.id}"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "/data/nhecker/apptainer/images/scanpy_1.11.4_coreinf_0.1.sif"

    input:
    tuple val(meta), path(h5ad)
    val(uns_key)
    val(species)

    output:
    tuple val(meta), path("*.h5ad"), emit: h5ad
    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    prefix = task.ext.prefix ?: "${meta.id}"
    template('enrich.py')
}
