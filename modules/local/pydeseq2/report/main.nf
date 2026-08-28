process PYDESEQ2_GENERATE_REPORT {
    tag "${meta.id}"
    label 'process_medium'

    container "${params.singularity_cache_dir
        ? params.singularity_cache_dir + '/scanpy_1.11.4_coreinf_0.3.sif'
        : '/nfsdata/apptainer/scanpy_1.11.4_coreinf_0.3.sif'}"

    input:
    tuple val(meta), path(h5ad)
    path(tsvs)
    path(manifest)
    path(qmd_template)

    output:
    tuple val(meta), path("*.html"), emit: html

    when:
    task.ext.when == null || task.ext.when

    script:
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    export NUMBA_CACHE_DIR=./tmp/numba
    export MPLCONFIGDIR=./tmp/matplotlib
    export XDG_CACHE_HOME=./tmp/matplotlib/cache

    export HOME=\$(readlink -f .)

    quarto render ${qmd_template} \\
        --output ${prefix}.html \\
        -P PSEUDOBULK_FILE:${h5ad} \\
        -P MANIFEST:${manifest} \\
        --to html
    """

    stub:
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch "${prefix}.html"
    """
}
