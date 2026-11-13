process MARKERS_COSG {
    tag "${meta.id}"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    // Built from containers/scanpy-ukdri-suite/Dockerfile.
    container 'docker.io/nhecker/scanpy-ukdri-suite:1.11.4-1'

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
        scanpy: \$(python3 -c 'import scanpy; print(scanpy.__version__)')
        cosg: \$(python3 -c 'from importlib.metadata import version; print(version("cosg"))')
        anndata: \$(python3 -c 'import anndata; print(anndata.__version__)')
        numpy: \$(python3 -c 'import numpy; print(numpy.__version__)')
        pandas: \$(python3 -c 'import pandas; print(pandas.__version__)')
    END_VERSIONS
    """
}
