process MARKERS_PRESTO {
    tag "${meta.id}"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    // Public base image, no derived image needed: this template only uses scanpy, anndata,
    // numpy, scipy, pandas, yaml and threadpoolctl, all of which the base already ships. Under
    // singularity/apptainer Nextflow converts it to a .sif in NXF_SINGULARITY_CACHEDIR.
    container 'docker.io/gcfntnu/scanpy:1.11.4'

    input:
    tuple val(meta), path(h5ad)
    val(obs_key)

    output:
    tuple val(meta), path("*.h5ad"), emit: h5ad, optional: true
    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    prefix = task.ext.prefix ?: "${meta.id}_markers_presto"
    template('presto.py')

    stub:
    prefix = task.ext.prefix ?: "${meta.id}_markers_presto"
    """
    touch "${prefix}.h5ad"

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python3 -c 'import platform; print(platform.python_version())')
        scanpy: \$(python3 -c 'import scanpy; print(scanpy.__version__)')
        anndata: \$(python3 -c 'import anndata; print(anndata.__version__)')
        numpy: \$(python3 -c 'import numpy; print(numpy.__version__)')
        scipy: \$(python3 -c 'import scipy; print(scipy.__version__)')
        pandas: \$(python3 -c 'import pandas; print(pandas.__version__)')
    END_VERSIONS
    """
}
