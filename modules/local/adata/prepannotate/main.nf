process ADATA_PREPANNOTATE {
    tag "${meta.id}"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    // Public base image: anndata, numpy, pandas and scipy are all this step needs.
    container 'docker.io/gcfntnu/scanpy:1.11.4'

    input:
    tuple val(meta), path(h5ad)
    val batch_col  // obs column with the query batches; one batch for all cells if missing
    val symbol_col // var column with the gene symbols, or 'index' for the var_names
    val subset_col    // obs column to select cells by; '' for all cells
    val subset_values // comma-separated values of subset_col to keep; '' for all cells
    val prediction_cols // comma-separated prediction columns of this run, e.g. 'celltypist:Adult_Human_Brain,scanvi:my_ref'

    output:
    tuple val(meta), path("${prefix}.h5ad"), emit: h5ad
    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    prefix = task.ext.prefix ?: "${meta.id}_annotate_input"

    if ("${h5ad}" == "${prefix}.h5ad") {
        error("Input and output names are the same, set prefix in module configuration to disambiguate!")
    }
    template('prepannotate.py')

    stub:
    prefix = task.ext.prefix ?: "${meta.id}_annotate_input"
    """
    touch ${prefix}.h5ad

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python3 -c 'import platform; print(platform.python_version())')
        anndata: \$(python3 -c 'import anndata; print(anndata.__version__)')
        numpy: \$(python3 -c 'import numpy; print(numpy.__version__)')
        pandas: \$(python3 -c 'import pandas; print(pandas.__version__)')
        scipy: \$(python3 -c 'import scipy; print(scipy.__version__)')
    END_VERSIONS
    """
}
