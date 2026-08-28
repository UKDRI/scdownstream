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
    val(sample_col)
    val(group_cols)
    val(ntop)
    path(versions_yaml)

    output:
    tuple val(meta), path("*.html"), emit: html
    path "versions.yml", emit: versions

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
        -P SAMPLE_COLUMN:"${sample_col}" \\
        -P GROUP_COLUMNS:"${group_cols}" \\
        -P NTOP:${ntop} \\
        -P VERSIONS_FILE:${versions_yaml} \\
        --to html

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        quarto: \$(quarto --version)
        python: \$(python3 -c 'import platform; print(platform.python_version())')
        anndata: \$(python3 -c 'import anndata; print(anndata.__version__)')
        matplotlib: \$(python3 -c 'import matplotlib; print(matplotlib.__version__)')
        numpy: \$(python3 -c 'import numpy; print(numpy.__version__)')
        pandas: \$(python3 -c 'import pandas; print(pandas.__version__)')
        ipython: \$(python3 -c 'import IPython; print(IPython.__version__)')
    END_VERSIONS
    """

    stub:
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch "${prefix}.html"
    touch "versions.yml"
    """
}
