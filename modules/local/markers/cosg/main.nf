process MARKERS_COSG {
    tag "${meta.id}"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    // BioContainer of bioconda's cosg: has anndata, numpy, pandas, scipy, scikit-learn and numba,
    // but no scanpy and no PyYAML, so the template uses neither.
    container "${workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container
        ? 'https://depot.galaxyproject.org/singularity/cosg:1.2.1--pyhdfd78af_0'
        : 'quay.io/biocontainers/cosg:1.2.1--pyhdfd78af_0'}"

    input:
    tuple val(meta), path(h5ad)
    val(obs_key)
    val(mu)
    val(rank_by)
    val(presto_ran)

    output:
    tuple val(meta), path("${prefix}.h5ad"), emit: h5ad
    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    prefix = task.ext.prefix ?: "${meta.id}_markers_cosg"
    template('cosg.py')

    stub:
    prefix = task.ext.prefix ?: "${meta.id}_markers_cosg"
    """
    touch "${prefix}.h5ad"

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python3 -c 'import platform; print(platform.python_version())')
        cosg: \$(python3 -c 'from importlib.metadata import version; print(version("cosg"))')
        anndata: \$(python3 -c 'import anndata; print(anndata.__version__)')
        numpy: \$(python3 -c 'import numpy; print(numpy.__version__)')
        pandas: \$(python3 -c 'import pandas; print(pandas.__version__)')
    END_VERSIONS
    """
}
