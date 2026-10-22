include { ADATA_TORDS           } from '../../modules/local/adata/tords'
include { ADATA_PUBLISH         } from '../../modules/local/adata/publish'
include { SCANPY_EXPORT_MARKERS } from '../../modules/local/scanpy/export_markers'


workflow FINALIZE_H5AD_SIMPLE {
    take:
        ch_h5ad

    main:
    ch_versions = Channel.empty()

    ADATA_TORDS(ch_h5ad)
    ch_versions = ch_versions.mix(ADATA_TORDS.out.versions)

    ADATA_PUBLISH(ch_h5ad)

    emit:
    ch_h5ad
    versions = ch_versions // channel: [ versions.yml ]
}


workflow FINALIZE_H5AD {
    take:
        ch_h5ad

    main:
    ch_versions = Channel.empty()

    ADATA_TORDS(ch_h5ad)
    ch_versions = ch_versions.mix(ADATA_TORDS.out.versions)

    ADATA_PUBLISH(ch_h5ad)
    SCANPY_EXPORT_MARKERS(ch_h5ad, params.name ? params.name : "markers", params.markers_uns_key, params.markers_thr_adj_pvalue, params.markers_n_top, params.markers_pct_nz, params.markers_min_logfc, params.markers_sep)

    emit:
    ch_h5ad
    json = SCANPY_EXPORT_MARKERS.out.json
    versions = ch_versions // channel: [ versions.yml ]
}
