process SCANPY_SCRUBLET {
    tag "${meta.id}"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container
        ? 'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/cc/cc499354bb89be6429cfc8d8416a340c2327a72732fb6ae33c346feef17b5ffb/data'
        : 'community.wave.seqera.io/library/python_pyyaml_scanpy_scikit-image:750e7b74b6d036e4'}"

    input:
    tuple val(meta), path(h5ad)
    val batch_col

    output:
    tuple val(meta), path("*.h5ad"), emit: h5ad
    tuple val(meta), path("*.pkl") , emit: predictions
    path "versions.yml"            , emit: versions
    path ("*_mqc.json"), emit: multiqc_files

    when:
    task.ext.when == null || task.ext.when

    script:
    prefix = task.ext.prefix ?: "${meta.id}_scrublet"
    if ("${h5ad}" == "${prefix}.h5ad") {
        error("Input and output names are the same, use \"task.ext.prefix\" to disambiguate!")
    }
    section_name = task.ext.section_name ?: "Scrublet"
    description = task.ext.description ?: "Doublet detection using scrublet. Scrublet simulates doublets by combining gene expression profiles of pairs of cells and trains a KNN classifier for predicting doublets."
    template('scrublet.py')

    stub:
    prefix = task.ext.prefix ?: "${meta.id}_scrublet"
    if ("${h5ad}" == "${prefix}.h5ad") {
        error("Input and output names are the same, use \"task.ext.prefix\" to disambiguate!")
    }
    """
    # These are needed to prevent errors during import of scanpy
    # when using singularity/apptainer
    export MPLCONFIGDIR=./tmp/mpl
    export NUMBA_CACHE_DIR=./tmp/numba

    touch ${prefix}.h5ad
    touch ${prefix}.pkl
    touch ${prefix}_scrublet_distributions.png
    touch ${prefix}_mqc.json

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python3 -c 'import platform; print(platform.python_version())')
        scanpy: \$(python3 -c 'import scanpy; print(scanpy.__version__)')
    END_VERSIONS
    """
}
