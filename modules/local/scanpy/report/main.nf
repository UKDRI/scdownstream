process SCANPY_GENERATE_REPORT {
    tag "${meta.id}"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    // Built from containers/scdownstream-tools/Dockerfile.
    container 'docker.io/nhecker/scdownstream-tools:1.11.4-1'

    input:
    tuple val(meta), path(h5ad)
    path(ipynb_template)
    val(clustering_name)
    val(ntop)
    val(umap_for_plots)
    val(umap_color_by)
    val(umap_color_by_embeddings)
    val(max_celltype_labels)

    output:
    tuple val(meta), path("*.html"), emit: html
    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    prefix = task.ext.prefix ?: "${meta.id}"
    // Only pass the colouring parameters when set; the report defaults cover the empty case
    color_by_args = [
        umap_color_by ? "-P \"UMAP_COLOR_BY:${umap_color_by}\"" : '',
        umap_color_by_embeddings ? "-P \"UMAP_COLOR_BY_EMBEDDINGS:${umap_color_by_embeddings}\"" : '',
    ].findAll().join(' ')
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
        -P "UMAP_FOR_PLOTS:${umap_for_plots}" \\
        -P MAX_CELLTYPE_LABELS:${max_celltype_labels} \\
        ${color_by_args} \\
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

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        quarto: \$(quarto --version)
        python: \$(python3 -c 'import platform; print(platform.python_version())')
        scanpy: \$(python3 -c 'import scanpy; print(scanpy.__version__)')
        anndata: \$(python3 -c 'import anndata; print(anndata.__version__)')
        liana: \$(python3 -c 'import liana; print(liana.__version__)')
    END_VERSIONS
    """
}


process SCANPY_GENERATE_REPORT_QC {
    tag "${meta.id}"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    // Built from containers/scdownstream-tools/Dockerfile.
    container 'docker.io/nhecker/scdownstream-tools:1.11.4-1'

    input:
    tuple val(meta), path(h5ad)
    path(ipynb_template)
    val(umap_for_plots)
    val(umap_color_by)
    val(umap_color_by_embeddings)
    val(max_celltype_labels)

    output:
    tuple val(meta), path("*.html"), emit: html
    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    prefix = task.ext.prefix ?: "${meta.id}"
    // Only pass the colouring parameters when set; the report defaults cover the empty case
    color_by_args = [
        umap_color_by ? "-P \"UMAP_COLOR_BY:${umap_color_by}\"" : '',
        umap_color_by_embeddings ? "-P \"UMAP_COLOR_BY_EMBEDDINGS:${umap_color_by_embeddings}\"" : '',
    ].findAll().join(' ')
    """
    export NUMBA_CACHE_DIR=./tmp/numba
    export MPLCONFIGDIR=./tmp/matplotlib
    export XDG_CACHE_HOME=./tmp/matplotlib/cache

    export HOME=\$(readlink -f .)

    quarto render ${ipynb_template} \\
        --output ${prefix}.html \\
        -P FILE:${h5ad} \\
        -P "UMAP_FOR_PLOTS:${umap_for_plots}" \\
        -P MAX_CELLTYPE_LABELS:${max_celltype_labels} \\
        ${color_by_args} \\
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

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        quarto: \$(quarto --version)
        python: \$(python3 -c 'import platform; print(platform.python_version())')
        scanpy: \$(python3 -c 'import scanpy; print(scanpy.__version__)')
        anndata: \$(python3 -c 'import anndata; print(anndata.__version__)')
    END_VERSIONS
    """
}


process SCANPY_GENERATE_REPORT_TRAINING {
    tag "${meta.id}"
    label 'process_low'

    conda "${moduleDir}/environment.yml"
    // Built from containers/scdownstream-tools/Dockerfile.
    container 'docker.io/nhecker/scdownstream-tools:1.11.4-1'

    input:
    tuple val(meta), path(label_counts)
    path(evaluation_files, stageAs: 'evaluation/*') // <method>_metrics/_per_label/_confusion.tsv; may be empty
    path(qmd_template)
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
        -P LABEL_COUNTS:${label_counts} \\
        -P EVALUATION_DIR:evaluation \\
        -P "MODEL_NAME:${meta.id}" \\
        -P VERSIONS_FILE:${versions_yaml} \\
        --to html

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        quarto: \$(quarto --version)
        python: \$(python3 -c 'import platform; print(platform.python_version())')
        matplotlib: \$(python3 -c 'import matplotlib; print(matplotlib.__version__)')
        numpy: \$(python3 -c 'import numpy; print(numpy.__version__)')
        pandas: \$(python3 -c 'import pandas; print(pandas.__version__)')
    END_VERSIONS
    """

    stub:
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch "${prefix}.html"

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        quarto: \$(quarto --version)
        python: \$(python3 -c 'import platform; print(platform.python_version())')
        matplotlib: \$(python3 -c 'import matplotlib; print(matplotlib.__version__)')
        numpy: \$(python3 -c 'import numpy; print(numpy.__version__)')
        pandas: \$(python3 -c 'import pandas; print(pandas.__version__)')
    END_VERSIONS
    """
}


process SCANPY_GENERATE_REPORT_ANNOTATION {
    tag "${meta.id}"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    // Built from containers/scdownstream-tools/Dockerfile.
    container 'docker.io/nhecker/scdownstream-tools:1.11.4-1'

    input:
    tuple val(meta), path(h5ad)
    path(qmd_template)
    val(prediction_cols) // comma-separated prediction columns added by this run, e.g. celltypist:<model>,scanvi:<model>
    val(subset)          // description of the cell subset, e.g. "leiden: 3, 7"; '' for all cells
    val(umap_for_plots)
    val(max_celltype_labels)
    path(versions_yaml)

    output:
    tuple val(meta), path("*.html"), emit: html
    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    prefix = task.ext.prefix ?: "${meta.id}"
    // quarto reads a -P value as YAML: the values are passed as single-quoted YAML strings, so the
    // colons in the column names and the subset ("leiden: 3, 7") stay text
    prediction_cols_yaml = "'" + (prediction_cols ?: '').toString().replace("'", "''") + "'"
    subset_yaml = "'" + (subset ?: '').toString().replace("'", "''") + "'"
    """
    export NUMBA_CACHE_DIR=./tmp/numba
    export MPLCONFIGDIR=./tmp/matplotlib
    export XDG_CACHE_HOME=./tmp/matplotlib/cache

    export HOME=\$(readlink -f .)

    quarto render ${qmd_template} \\
        --output ${prefix}.html \\
        -P FILE:${h5ad} \\
        -P "PREDICTION_COLS:${prediction_cols_yaml}" \\
        -P "SUBSET:${subset_yaml}" \\
        -P "UMAP_FOR_PLOTS:${umap_for_plots}" \\
        -P MAX_CELLTYPE_LABELS:${max_celltype_labels} \\
        -P VERSIONS_FILE:${versions_yaml} \\
        --to html

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        quarto: \$(quarto --version)
        python: \$(python3 -c 'import platform; print(platform.python_version())')
        scanpy: \$(python3 -c 'import scanpy; print(scanpy.__version__)')
        anndata: \$(python3 -c 'import anndata; print(anndata.__version__)')
        matplotlib: \$(python3 -c 'import matplotlib; print(matplotlib.__version__)')
        pandas: \$(python3 -c 'import pandas; print(pandas.__version__)')
    END_VERSIONS
    """

    stub:
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch "${prefix}.html"

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        quarto: \$(quarto --version)
        python: \$(python3 -c 'import platform; print(platform.python_version())')
        scanpy: \$(python3 -c 'import scanpy; print(scanpy.__version__)')
        anndata: \$(python3 -c 'import anndata; print(anndata.__version__)')
        matplotlib: \$(python3 -c 'import matplotlib; print(matplotlib.__version__)')
        pandas: \$(python3 -c 'import pandas; print(pandas.__version__)')
    END_VERSIONS
    """
}
