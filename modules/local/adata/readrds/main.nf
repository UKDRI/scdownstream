process ADATA_READRDS {
    tag "$meta.id"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    // The base image ships rpy2, anndata2ri (git main, >= 2.0), r-seurat and bioconductor-singlecellexperiment
    container 'docker.io/gcfntnu/scanpy:1.11.4'

    input:
    tuple val(meta), path(rds)

    output:
    tuple val(meta), path("*.h5ad"), emit: h5ad
    path "versions.yml"            , emit: versions

    script:
    prefix = task.ext.prefix ?: "${meta.id}"
    template 'readrds.py'

    stub:
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch "${prefix}.h5ad"
    touch "versions.yml"
    """
}
