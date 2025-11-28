/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT MODULES / SUBWORKFLOWS / FUNCTIONS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { LOAD_H5AD                            } from '../subworkflows/local/load_h5ad'
include { QUALITY_CONTROL                      } from '../subworkflows/local/quality_control'
include { UNIFY                                } from '../subworkflows/local/unify'
include { CELLTYPE_ASSIGNMENT                  } from '../subworkflows/local/celltype_assignment'
include { ADATA_EXTEND as FINALIZE_QC_ANNDATAS } from '../modules/local/adata/extend'
include { COMBINE                              } from '../subworkflows/local/combine'
include { ADATA_SPLITEMBEDDINGS                } from '../modules/local/adata/splitembeddings'
include { CLUSTER                              } from '../subworkflows/local/cluster'
include { PSEUDOBULKING                        } from '../subworkflows/local/pseudobulking'
include { PER_GROUP                            } from '../subworkflows/local/per_group'
include { FINALIZE                             } from '../subworkflows/local/finalize'
include { MULTIQC                              } from '../modules/nf-core/multiqc/main'
include { paramsSummaryMap                     } from 'plugin/nf-schema'
include { paramsSummaryMultiqc                 } from '../subworkflows/nf-core/utils_nfcore_pipeline'
include { softwareVersionsToYAML               } from '../subworkflows/nf-core/utils_nfcore_pipeline'
include { methodsDescriptionText               } from '../subworkflows/local/utils_nfcore_scdownstream_pipeline'
include { SCANPY_HVGS                          } from '../modules/local/scanpy/hvgs'
include { SCANPY_PCA                           } from '../modules/local/scanpy/pca'
include { SCANPY_NEIGHBORS                     } from '../modules/local/scanpy/neighbors'
include { SCANPY_UMAP                          } from '../modules/local/scanpy/umap'
include { SCANPY_LOG_NORMALIZE                 } from '../modules/local/scanpy/normalization'
include { SCANPY_LEIDEN                        } from '../modules/local/scanpy/leiden'
include { SCANPY_RANKGENESGROUPS               } from '../modules/local/scanpy/rankgenesgroups'
include { SCANPY_ENRICH                        } from '../modules/local/scanpy/enrich'
include { SCANPY_GENERATE_REPORT               } from '../modules/local/scanpy/report'
include { SCANPY_REPORT_TO_HTML                } from '../modules/local/scanpy/report'
include { LIANA_RANKAGGREGATE                  } from '../modules/local/liana/rankaggregate'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN MAIN WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow SCDOWNSTREAM {
    take:
    ch_samplesheet // channel: samplesheet read in from --input
    ch_base        // channel: [ val(meta), path(h5ad) ]

    main:

    ch_versions = Channel.empty()
    ch_integrations = Channel.empty()
    ch_obs = Channel.empty()
    ch_var = Channel.empty()
    ch_obsm = Channel.empty()
    ch_obsp = Channel.empty()
    ch_uns = Channel.empty()
    ch_layers = Channel.empty()
    ch_multiqc_files = Channel.empty()
    ch_h5ad = Channel.empty()
    ch_combined = Channel.empty()


    ch_obs_per_sample = Channel.empty()
    ch_var_per_sample = Channel.empty()
    ch_obsm_per_sample = Channel.empty()
    ch_obsp_per_sample = Channel.empty()
    ch_uns_per_sample = Channel.empty()
    ch_layers_per_sample = Channel.empty()

    //
    // Load/Convert input to h5ad
    //
    LOAD_H5AD(ch_samplesheet)
    ch_h5ad = LOAD_H5AD.out.h5ad
    ch_versions = ch_versions.mix(LOAD_H5AD.out.versions)

    //
    // Quality control per sample
    //
    QUALITY_CONTROL(
        ch_h5ad,
        params.ambient_correction,
        !params.doublet_detection || params.doublet_detection == 'none' ? [] : params.doublet_detection.split(',').collect { it -> it.trim().toLowerCase() },
        params.doublet_detection_threshold,
        params.mito_genes,
    )
    ch_versions = ch_versions.mix(QUALITY_CONTROL.out.versions)
    ch_multiqc_files = ch_multiqc_files.mix(QUALITY_CONTROL.out.multiqc_files)
    ch_h5ad = QUALITY_CONTROL.out.h5ad

    //
    // Perform automated celltype assignment
    //
    CELLTYPE_ASSIGNMENT(ch_h5ad.map { meta, h5ad -> [meta, h5ad, meta.symbol_col] })
    ch_versions = ch_versions.mix(CELLTYPE_ASSIGNMENT.out.versions)
    ch_obs_per_sample = ch_obs_per_sample.mix(CELLTYPE_ASSIGNMENT.out.obs)

    FINALIZE_QC_ANNDATAS(
        ch_h5ad.join(ch_obs_per_sample.groupTuple(), remainder: true).join(ch_var_per_sample.groupTuple(), remainder: true).join(ch_obsm_per_sample.groupTuple(), remainder: true).join(ch_obsp_per_sample.groupTuple(), remainder: true).join(ch_uns_per_sample.groupTuple(), remainder: true).join(ch_layers_per_sample.groupTuple(), remainder: true).map { meta, h5ad, obs, var, obsm, obsp, uns, layers ->
            [meta, h5ad, obs ?: [], var ?: [], obsm ?: [], obsp ?: [], uns ?: [], layers ?: []]
        }
    )
    ch_h5ad = FINALIZE_QC_ANNDATAS.out.h5ad
    ch_versions = ch_versions.mix(FINALIZE_QC_ANNDATAS.out.versions)

    if (!params.qc_only) {
        //
        // Unify samples to make them compatible for integration
        //
        UNIFY(ch_h5ad)
        ch_versions = ch_versions.mix(UNIFY.out.versions)
        ch_multiqc_files = ch_multiqc_files.mix(UNIFY.out.multiqc_files)
        ch_h5ad = UNIFY.out.h5ad

        //
        // Combine samples and perform integration
        //
        COMBINE(ch_h5ad, ch_base)
        ch_versions = ch_versions.mix(COMBINE.out.versions)
        ch_obs = ch_obs.mix(COMBINE.out.obs)
        ch_var = ch_var.mix(COMBINE.out.var)
        ch_obsm = ch_obsm.mix(COMBINE.out.obsm)
        ch_integrations = ch_integrations.mix(COMBINE.out.integrations)
        //ch_finalization_base = COMBINE.out.h5ad
        ch_combined = COMBINE.out.h5ad

        ch_label_grouping = COMBINE.out.h5ad_inner
        grouping_col = "label"
    }

    //
    // Compute embeddings
    //
    if (!params.qc_only) {

        ch_h5ad = ch_combined
        SCANPY_LOG_NORMALIZE(ch_h5ad)
        ch_h5ad = SCANPY_LOG_NORMALIZE.out.h5ad
        SCANPY_HVGS(ch_h5ad, params.n_hvgs, false)
        ch_h5ad = SCANPY_HVGS.out.h5ad
        SCANPY_PCA(ch_h5ad)
        ch_h5ad = SCANPY_PCA.out.h5ad
        SCANPY_NEIGHBORS(ch_h5ad, 'X_pca', 'neighbors_pca' )
        ch_h5ad = SCANPY_NEIGHBORS.out.h5ad
        SCANPY_UMAP(ch_h5ad, 'neighbors_pca', 'X_umap')
        ch_h5ad = SCANPY_UMAP.out.h5ad

    }


    //
    // Perform clustering and related analysis
    //
    if (!params.qc_only) {
        SCANPY_LEIDEN(ch_h5ad, params.clustering_resolution, "leiden", params.cluster_neighbors, false)
        ch_h5ad = SCANPY_LEIDEN.out.h5ad
        SCANPY_RANKGENESGROUPS(ch_h5ad, "leiden")
        ch_h5ad = SCANPY_RANKGENESGROUPS.out.h5ad
        SCANPY_ENRICH(ch_h5ad, "rank_genes_groups", params.species, params.enrich_min_in_group_fraction, params.enrich_min_fold_change, params.enrich_max_out_group_fraction)
        ch_h5ad = SCANPY_ENRICH.out.h5ad
        LIANA_RANKAGGREGATE(ch_h5ad, params.species)
        ch_h5ad = LIANA_RANKAGGREGATE.out.h5ad
    }

    //
    // Summary report
    //
    if (!params.qc_only) {
        SCANPY_GENERATE_REPORT(ch_h5ad, "${projectDir}/modules/local/scanpy/report/templates/scdownstream_report.ipynb")
        SCANPY_REPORT_TO_HTML(SCANPY_GENERATE_REPORT.out.ipynp)
    }

    //
    // Perform clustering and per-cluster analysis
    //
    // if (!params.qc_only) {
    //    CLUSTER(
    //        ch_integrations,
    //        params.cluster_per_label,
    //        params.cluster_global,
    //        params.input ? "label" : params.base_label_col,
    //        params.clustering_resolutions.split(','),
    //        "batch",
    //        "X_emb",
    //    )
    //    ch_versions = ch_versions.mix(CLUSTER.out.versions)
    //    ch_obs = ch_obs.mix(CLUSTER.out.obs)
    //    ch_obsm = ch_obsm.mix(CLUSTER.out.obsm)
    //    ch_multiqc_files = ch_multiqc_files.mix(CLUSTER.out.multiqc_files)

    //    if (params.pseudobulk) {
    //        PSEUDOBULKING(
    //            CLUSTER.out.h5ad_clustering,
    //            params.pseudobulk_groupby_labels.split(','),
    //            params.pseudobulk_min_num_cells,
    //            "X",
    //        )
    //        ch_versions = ch_versions.mix(PSEUDOBULKING.out.versions)
    //    }

    //    PER_GROUP(
    //        CLUSTER.out.h5ad_clustering.map { meta, h5ad -> [meta + [obs_key: "${meta.id}_leiden"], h5ad] },
    //        CLUSTER.out.h5ad_neighbors.map { meta, h5ad -> [meta + [obs_key: grouping_col], h5ad] },
    //        ch_label_grouping.map { meta, h5ad -> [meta + [obs_key: grouping_col], h5ad] },
    //    )
    //    ch_versions = ch_versions.mix(PER_GROUP.out.versions)
    //    ch_uns = ch_uns.mix(PER_GROUP.out.uns)
    //    ch_multiqc_files = ch_multiqc_files.mix(PER_GROUP.out.multiqc_files)

    //    FINALIZE(ch_finalization_base, ch_obs, ch_var, ch_obsm, ch_obsp, ch_uns, ch_layers)
    //    ch_versions = ch_versions.mix(FINALIZE.out.versions)
    //}

    //
    // Collate and save software versions
    //
    softwareVersionsToYAML(ch_versions)
        .collectFile(
            storeDir: "${params.outdir}/pipeline_info",
            name: 'nf_core_' + 'scdownstream_software_' + 'mqc_' + 'versions.yml',
            sort: true,
            newLine: true,
        )
        .set { ch_collated_versions }

    //
    // MODULE: MultiQC
    //
    ch_multiqc_config = Channel.fromPath(
        "${projectDir}/assets/multiqc_config.yml",
        checkIfExists: true
    )
    ch_multiqc_custom_config = params.multiqc_config
        ? Channel.fromPath(params.multiqc_config, checkIfExists: true)
        : Channel.empty()
    ch_multiqc_logo = params.multiqc_logo
        ? Channel.fromPath(params.multiqc_logo, checkIfExists: true)
        : Channel.empty()

    summary_params = paramsSummaryMap(
        workflow,
        parameters_schema: "nextflow_schema.json"
    )
    ch_workflow_summary = Channel.value(paramsSummaryMultiqc(summary_params))
    ch_multiqc_files = ch_multiqc_files.mix(
        ch_workflow_summary.collectFile(name: 'workflow_summary_mqc.yaml')
    )
    ch_multiqc_custom_methods_description = params.multiqc_methods_description
        ? file(params.multiqc_methods_description, checkIfExists: true)
        : file("${projectDir}/assets/methods_description_template.yml", checkIfExists: true)
    ch_methods_description = Channel.value(
        methodsDescriptionText(ch_multiqc_custom_methods_description)
    )

    ch_multiqc_files = ch_multiqc_files.mix(ch_collated_versions)
    ch_multiqc_files = ch_multiqc_files.mix(
        ch_methods_description.collectFile(
            name: 'methods_description_mqc.yaml',
            sort: true,
        )
    )

    MULTIQC(
        ch_multiqc_files.collect(),
        ch_multiqc_config.toList(),
        ch_multiqc_custom_config.toList(),
        ch_multiqc_logo.toList(),
        [],
        [],
    )

    emit:
    multiqc_report = MULTIQC.out.report.toList() // channel: /path/to/multiqc_report.html
    versions       = ch_versions // channel: [ path(versions.yml) ]
}
