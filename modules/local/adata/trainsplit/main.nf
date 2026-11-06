process ADATA_TRAINSPLIT {
    tag "${meta.id}"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    // Public base image: anndata, pandas and scikit-learn are all this step needs.
    container 'docker.io/gcfntnu/scanpy:1.11.4'

    input:
    tuple val(meta), path(h5ad)

    output:
    tuple val(meta), path("${prefix}_all.h5ad"), emit: all
    tuple val(meta), path("${prefix}_train.h5ad"), emit: train, optional: true
    tuple val(meta), path("${prefix}_test.h5ad"), emit: test, optional: true
    tuple val(meta), path("${prefix}_label_counts.tsv"), emit: label_counts
    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    prefix = task.ext.prefix ?: "${meta.id}"
    label_col = task.ext.label_col ?: ''
    unknown_label = task.ext.unknown_label ?: 'unknown'
    batch_col = task.ext.batch_col ?: 'batch'
    counts_layer = task.ext.counts_layer ?: 'X'
    symbol_col = task.ext.symbol_col ?: 'index'
    min_cells_per_label = task.ext.min_cells_per_label != null ? task.ext.min_cells_per_label : 50
    test_fraction = task.ext.test_fraction ?: 0.25
    split_by = task.ext.split_by ?: ''
    evaluate = task.ext.evaluate != null ? task.ext.evaluate : true
    // null or 0: no cap
    max_cells_per_label = task.ext.max_cells_per_label ?: 0

    if ("${h5ad}" in ["${prefix}_all.h5ad", "${prefix}_train.h5ad", "${prefix}_test.h5ad"]) {
        error("Input and output names are the same, set prefix in module configuration to disambiguate!")
    }
    template('trainsplit.py')

    stub:
    prefix = task.ext.prefix ?: "${meta.id}"
    evaluate = task.ext.evaluate != null ? task.ext.evaluate : true
    """
    touch ${prefix}_all.h5ad
    touch ${prefix}_label_counts.tsv

    if [ "${evaluate}" = "true" ]; then
        touch ${prefix}_train.h5ad
        touch ${prefix}_test.h5ad
    fi

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python3 -c 'import platform; print(platform.python_version())')
        anndata: \$(python3 -c 'import anndata; print(anndata.__version__)')
        numpy: \$(python3 -c 'import numpy; print(numpy.__version__)')
        pandas: \$(python3 -c 'import pandas; print(pandas.__version__)')
        scipy: \$(python3 -c 'import scipy; print(scipy.__version__)')
        scikit-learn: \$(python3 -c 'import sklearn; print(sklearn.__version__)')
    END_VERSIONS
    """
}
