include { CELLTYPES_CELLTYPIST } from '../../../modules/local/celltypes/celltypist'
include { ADATA_ADD_OBS_OBSM   } from '../../../modules/local/adata/add_obs_obsm'

workflow CELLTYPE_ASSIGNMENT {
    take:
    ch_h5ad // channel: [ meta, h5ad, symbol_col ]

    main:
    ch_versions = Channel.empty()
    ch_obs = Channel.empty()

    if (params.celltypist_model) {
        celltypist_entries = params.celltypist_model.split(',').collect { entry -> entry.trim() }.findAll { entry -> entry }
        // Entries ending in .pkl or holding a path are local model files (e.g. from --stage
        // train_models); everything else is the name of a built-in model to download
        celltypist_names = celltypist_entries.findAll { entry -> !(entry.endsWith('.pkl') || entry.contains('/')) }
        celltypist_files = celltypist_entries
            .findAll { entry -> entry.endsWith('.pkl') || entry.contains('/') }
            .collect { entry -> file(entry, checkIfExists: true) }

        CELLTYPES_CELLTYPIST(ch_h5ad, celltypist_names, celltypist_files)

        ADATA_ADD_OBS_OBSM(ch_h5ad.join(CELLTYPES_CELLTYPIST.out.h5ad).map{ meta, h5ad1, _symbol_col, h5ad2 -> [meta, h5ad1, h5ad2]})
        ch_h5ad = ADATA_ADD_OBS_OBSM.out.h5ad
        ch_obs = ch_obs.mix(CELLTYPES_CELLTYPIST.out.obs)
        ch_versions = ch_versions.mix(CELLTYPES_CELLTYPIST.out.versions)
    }

    emit:
    h5ad     = ch_h5ad     // channel: [meta, h5ad]
    obs      = ch_obs      // channel: [ meta, pkl ]
    versions = ch_versions // channel: [ versions.yml ]
}
