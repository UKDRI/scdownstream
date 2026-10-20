process PYDESEQ2_EXPORT_MARKERS {
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    // Public base image, no derived image needed: this template only uses pandas, numpy
    // and yaml, all of which the base already ships. Under singularity/apptainer Nextflow
    // converts it to a .sif in NXF_SINGULARITY_CACHEDIR.
    container 'docker.io/gcfntnu/scanpy:1.11.4'

    input:
    path(tsvs)
    path(manifest)
    val(project_name)
    val(group_col)
    val(thr_adj_pvalue)
    val(n_top)
    val(min_logfc)
    val(sep)

    output:
    path("${prefix}.json.gz"), emit: json
    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    prefix = task.ext.prefix ?: "differential_genes_markers"
    template('export_markers.py')

    stub:
    prefix = task.ext.prefix ?: "differential_genes_markers"
    """
    touch "${prefix}.json.gz"

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python3 -c 'import platform; print(platform.python_version())')
        pandas: \$(python3 -c 'import pandas; print(pandas.__version__)')
        numpy: \$(python3 -c 'import numpy; print(numpy.__version__)')
    END_VERSIONS
    """
}
