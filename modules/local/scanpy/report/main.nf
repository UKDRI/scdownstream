process SCANPY_GENERATE_REPORT {
    tag "${meta.id}"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "/nfsdata/apptainer/scanpy_1.11.4_coreinf_0.3.sif"

    input:
    tuple val(meta), path(h5ad)
    path(ipynb_template)
    val(clustering_name)
    val(ntop)
    
    output:
    tuple val(meta), path("*.html"), emit: html
    //path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    prefix = task.ext.prefix ?: "${meta.id}"
    //ipynb_template = "${moduleDir}/templates/scdownstream_report.ipynb"

    """
    export NUMBA_CACHE_DIR=./tmp/numba
    export MPLCONFIGDIR=./tmp/matplotlib
    export XDG_CACHE_HOME=./tmp/matplotlib/cache

    export HOME=`readlink -f .`

        quarto render ${ipynb_template} \
            --output ${prefix}.html \
            -P FILE:${h5ad} \
            -P CLUSTERING_NAME:${clustering_name} \
            -P NTOP:${ntop} \
            --to html
    """
}


process SCANPY_GENERATE_REPORT_QC {
    tag "${meta.id}"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "/nfsdata/apptainer/scanpy_1.11.4_coreinf_0.3.sif"

    input:
    tuple val(meta), path(h5ad)
    path(ipynb_template)
    
    output:
    tuple val(meta), path("*.html"), emit: html
    //path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    prefix = task.ext.prefix ?: "${meta.id}"
    //ipynb_template = "${moduleDir}/templates/scdownstream_report.ipynb"

    """
    export NUMBA_CACHE_DIR=./tmp/numba
    export MPLCONFIGDIR=./tmp/matplotlib
    export XDG_CACHE_HOME=./tmp/matplotlib/cache

    export HOME=`readlink -f .`

        quarto render ${ipynb_template} \
            --output ${prefix}.html \
            -P FILE:${h5ad} \
            --to html
    """
}

