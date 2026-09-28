process SCANPY_GENERATE_REPORT {
    tag "${meta.id}"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    // Built from modules/local/scanpy/report/Dockerfile.
    container 'docker.io/nhecker/scanpy-report:1.11.4-coreinf0.4'

    input:
    tuple val(meta), path(h5ad)
    path(ipynb_template)
    val(clustering_name)
    val(ntop)

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

    quarto render ${ipynb_template} \\
        --output ${prefix}.html \\
        -P FILE:${h5ad} \\
        -P CLUSTERING_NAME:${clustering_name} \\
        -P NTOP:${ntop} \\
        --to html

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        quarto: \$(quarto --version)
        python: \$(python3 -c 'import platform; print(platform.python_version())')
        scanpy: \$(python3 -c 'import scanpy; print(scanpy.__version__)')
        anndata: \$(python3 -c 'import anndata; print(anndata.__version__)')
        liana: \$(python3 -c 'import liana; print(liana.__version__)')
    END_VERSIONS
    """

    stub:
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch "${prefix}.html"
    touch "versions.yml"
    """
}


process SCANPY_GENERATE_REPORT_QC {
    tag "${meta.id}"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    // Built from modules/local/scanpy/report/Dockerfile.
    container 'docker.io/nhecker/scanpy-report:1.11.4-coreinf0.4'

    input:
    tuple val(meta), path(h5ad)
    path(ipynb_template)

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

    quarto render ${ipynb_template} \\
        --output ${prefix}.html \\
        -P FILE:${h5ad} \\
        --to html

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        quarto: \$(quarto --version)
        python: \$(python3 -c 'import platform; print(platform.python_version())')
        scanpy: \$(python3 -c 'import scanpy; print(scanpy.__version__)')
        anndata: \$(python3 -c 'import anndata; print(anndata.__version__)')
    END_VERSIONS
    """

    stub:
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch "${prefix}.html"
    touch "versions.yml"
    """
}
