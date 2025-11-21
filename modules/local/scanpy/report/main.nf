process SCANPY_GENERATE_REPORT {
    tag "${meta.id}"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "/data/nhecker/apptainer/images/scanpy_1.11.4_coreinf_0.1.sif"

    input:
    tuple val(meta), path(h5ad)
    path(ipynb_template)
    
    output:
    tuple val(meta), path("*.ipynb"), emit: ipynp
    //path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    prefix = task.ext.prefix ?: "${meta.id}"
    //ipynb_template = "${moduleDir}/templates/scdownstream_report.ipynb"

    """
    export NUMBA_CACHE_DIR=./tmp/numba
    export MPLCONFIGDIR=./tmp/matplotlib

    papermill ${ipynb_template} \
            --report-mode \
            ${prefix}.ipynb \
            -p FILE ${h5ad}
    """
}


process SCANPY_REPORT_TO_HTML {
    tag "${meta.id}"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "/data/nhecker/apptainer/images/scanpy_1.11.4_coreinf_0.1.sif"

    input:
    tuple val(meta), path(notebook)
    
    output:
    tuple val(meta), path("*.html"), emit: h5ad
    //path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    prefix = task.ext.prefix ?: "${meta.id}"

    """
    jupyter nbconvert ${notebook} --to html --no-input
    """
}
