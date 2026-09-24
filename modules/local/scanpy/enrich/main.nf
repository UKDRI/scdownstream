process SCANPY_ENRICH {
    tag "${meta.id}"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    // Cannot use the plain gcfntnu/scanpy:1.11.4 base: enrich.py calls sc.queries.enrich(),
    // which imports gprofiler-official, and the base does not ship it. The bare `except:` in
    // the template would swallow the ImportError and write 'empty set' for every group, so the
    // process would exit 0 with silently empty enrichment. Keep the derived image.
    // Built from modules/local/scanpy/report/Dockerfile.
    container 'docker.io/nhecker/scanpy-report:1.11.4-coreinf0.4'

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
