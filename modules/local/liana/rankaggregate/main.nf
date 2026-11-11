process LIANA_RANKAGGREGATE {
    tag "$meta.id"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    // The pipeline's report image (gcfntnu/scanpy:1.11.4 base, liana 1.6.1), built from
    // modules/local/scanpy/report/Dockerfile. The old BioContainer liana:1.4.0 had anndata 0.10,
    // which could not read the h5ad files written by the anndata 0.12 steps.
    container 'docker.io/nhecker/scanpy-report:1.11.4-coreinf0.4'

    input:
    tuple val(meta), path(h5ad)
    val(obs_key)
    val(species)
    path(hcop_dir)

    output:
    tuple val(meta), path("*.h5ad"), emit: h5ad, optional: true
    path("*.pkl")                  , emit: uns, optional: true
    path "versions.yml"            , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    prefix = task.ext.prefix ?: "${meta.id}"
    template 'rank_aggregate.py'

    stub:
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.h5ad
    touch ${prefix}.pkl
    touch versions.yml
    """
}
