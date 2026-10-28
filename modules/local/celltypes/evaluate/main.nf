process MODEL_EVALUATE {
    tag "${meta.id}"
    label 'process_low'

    conda "${moduleDir}/environment.yml"
    // Public base image: anndata, pandas and scikit-learn are all this step needs.
    container 'docker.io/gcfntnu/scanpy:1.11.4'

    input:
    tuple val(meta), path(h5ad)
    val pred_col
    val true_col

    output:
    tuple val(meta), path("${prefix}_metrics.tsv"), emit: metrics
    tuple val(meta), path("${prefix}_per_label.tsv"), emit: per_label
    tuple val(meta), path("${prefix}_confusion.tsv"), emit: confusion
    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    prefix = task.ext.prefix ?: "${meta.method}"
    template('evaluate.py')

    stub:
    prefix = task.ext.prefix ?: "${meta.method}"
    """
    touch ${prefix}_metrics.tsv
    touch ${prefix}_per_label.tsv
    touch ${prefix}_confusion.tsv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python3 -c 'import platform; print(platform.python_version())')
        anndata: \$(python3 -c 'import anndata; print(anndata.__version__)')
        numpy: \$(python3 -c 'import numpy; print(numpy.__version__)')
        pandas: \$(python3 -c 'import pandas; print(pandas.__version__)')
        scikit-learn: \$(python3 -c 'import sklearn; print(sklearn.__version__)')
    END_VERSIONS
    """
}
