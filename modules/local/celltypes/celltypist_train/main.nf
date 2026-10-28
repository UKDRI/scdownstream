process CELLTYPIST_TRAIN {
    tag "${meta.id}"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    // Same image as CELLTYPES_CELLTYPIST, so a model is trained and applied with the same versions
    container "${workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container
        ? 'oras://community.wave.seqera.io/library/celltypist_scanpy:89a98f51262cfff4'
        : 'community.wave.seqera.io/library/celltypist_scanpy:44b604b24dd4cf33'}"

    input:
    tuple val(meta), path(h5ad)

    output:
    tuple val(meta), path("${prefix}.pkl"), emit: model
    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    prefix = task.ext.prefix ?: "${meta.id}"
    feature_selection = task.ext.feature_selection != null ? task.ext.feature_selection : true
    use_sgd = task.ext.use_sgd ?: false
    mini_batch = task.ext.mini_batch ?: false
    batch_size = task.ext.batch_size ?: 1000
    template('celltypist_train.py')

    stub:
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.pkl

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python3 -c 'import platform; print(platform.python_version())')
        scanpy: \$(python3 -c 'import scanpy; print(scanpy.__version__)')
        celltypist: \$(python3 -c 'import celltypist; print(celltypist.__version__)')
        scikit-learn: \$(python3 -c 'import sklearn; print(sklearn.__version__)')
    END_VERSIONS
    """
}
