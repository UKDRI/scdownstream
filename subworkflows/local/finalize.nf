include { ADATA_EXTEND        } from '../../modules/local/adata/extend'
include { ADATA_TORDS         } from '../../modules/local/adata/tords'
include { ADATA_PREPCELLXGENE } from '../../modules/local/adata/prepcellxgene'
include { ADATA_PUBLISH         } from '../../modules/local/adata/publish'

workflow FINALIZE {
    take:
    ch_h5ad   // channel: [ merged, h5ad ]
    ch_obs    // channel: [ pkl ]
    ch_var    // channel: [ pkl ]
    ch_obsm   // channel: [ pkl ]
    ch_obsp
    ch_uns    // channel: [ pkl ]
    ch_layers

    main:
    ch_versions = Channel.empty()

    ADATA_EXTEND(ch_h5ad
        .combine(ch_obs.flatten().collect().ifEmpty([]).map{ it -> [it] })
        .combine(ch_var.flatten().collect().ifEmpty([]).map{ it -> [it] })
        .combine(ch_obsm.flatten().collect().ifEmpty([]).map{ it -> [it] })
        .combine(ch_obsp.flatten().collect().ifEmpty([]).map{ it -> [it] })
        .combine(ch_uns.flatten().collect().ifEmpty([]).map{ it -> [it] })
        .combine(ch_layers.flatten().collect().ifEmpty([]).map{ it -> [it] })
    )
    ch_versions = ch_versions.mix(ADATA_EXTEND.out.versions)

    ADATA_TORDS(ADATA_EXTEND.out.h5ad)
    ch_versions = ch_versions.mix(ADATA_TORDS.out.versions)

    if (params.prep_cellxgene) {
        ADATA_PREPCELLXGENE(ADATA_EXTEND.out.h5ad)
        ch_versions = ch_versions.mix(ADATA_PREPCELLXGENE.out.versions)
    }

    emit:
    versions = ch_versions // channel: [ versions.yml ]
}


workflow FINALIZE_H5AD {
    take:
        ch_h5ad

    main:
    ch_versions = Channel.empty()

    ADATA_TORDS(ch_h5ad)
    ch_versions = ch_versions.mix(ADATA_TORDS.out.versions)

    if (params.prep_cellxgene) {
        ADATA_PREPCELLXGENE(ADATA_EXTEND.out.h5ad)
        ch_versions = ch_versions.mix(ADATA_PREPCELLXGENE.out.versions)
    }

    ADATA_PUBLISH(ch_h5ad)

    emit:
    ch_h5ad
    versions = ch_versions // channel: [ versions.yml ]
}
