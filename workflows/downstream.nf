/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT MODULES / SUBWORKFLOWS / FUNCTIONS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/


include { MULTIQC                              } from '../modules/nf-core/multiqc/main'
include { paramsSummaryMap                     } from 'plugin/nf-schema'
include { paramsSummaryMultiqc                 } from '../subworkflows/nf-core/utils_nfcore_pipeline'
include { softwareVersionsToYAML               } from '../subworkflows/nf-core/utils_nfcore_pipeline'
include { methodsDescriptionText               } from '../subworkflows/local/utils'
include { MARKERS_PRESTO                       } from '../modules/local/markers/presto'
include { MARKERS_COSG                         } from '../modules/local/markers/cosg'
include { SCANPY_ENRICH                        } from '../modules/local/scanpy/enrich'
include { SCANPY_GENERATE_REPORT               } from '../modules/local/scanpy/report'
include { LIANA_RANKAGGREGATE                  } from '../modules/local/liana/rankaggregate'
include { FINALIZE_H5AD                        } from '../subworkflows/local/finalize'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN MAIN WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow DOWNSTREAM_ANALYSIS {
    take:
    ch_samplesheet // channel: samplesheet read in from --input
    ch_base        // channel: [ val(meta), path(h5ad) ]

    main:

    ch_versions = Channel.empty()

    //
    // Perform downstream analysis
    //
    clustering_name = params.selected_clustering == null ? 'leiden' :  params.selected_clustering

    ch_h5ad = ch_base

    // Marker gene methods from --markers_method (presto and/or cosg); the ranking defaults to
    // AUC when presto runs, else to the COSG score. validateMarkersParameters() checks both.
    def markers_methods = params.markers_method.toString().tokenize(',').collect { method -> method.trim().toLowerCase() }.findAll { method -> method }.unique()
    def markers_rank_by = params.markers_rank_by ?: ('presto' in markers_methods ? 'auc' : 'cosg')

    // MODULE: marker genes per group (presto Wilcoxon + AUC)
    if ('presto' in markers_methods) {
        MARKERS_PRESTO(ch_h5ad, clustering_name)
        ch_versions = ch_versions.mix(MARKERS_PRESTO.out.versions)
        ch_h5ad = MARKERS_PRESTO.out.h5ad
    }

    // MODULE: COSG score per group, added to the presto results if presto ran
    if ('cosg' in markers_methods) {
        MARKERS_COSG(
            ch_h5ad,
            clustering_name,
            params.markers_cosg_mu,
            markers_rank_by,
            'presto' in markers_methods
        )
        ch_versions = ch_versions.mix(MARKERS_COSG.out.versions)
        ch_h5ad = MARKERS_COSG.out.h5ad
    }
    SCANPY_ENRICH(ch_h5ad, "rank_genes_groups", params.species, params.enrich_min_in_group_fraction, params.enrich_min_fold_change, params.enrich_max_out_group_fraction)
    ch_h5ad = SCANPY_ENRICH.out.h5ad
    LIANA_RANKAGGREGATE(
        ch_h5ad,
        clustering_name,
        params.species,
        params.ortholog_hcop_directory ? file(params.ortholog_hcop_directory, checkIfExists: true) : []
    )
    ch_h5ad = LIANA_RANKAGGREGATE.out.h5ad

    // create output final output files
    FINALIZE_H5AD(ch_h5ad)

    ch_versions = ch_versions.mix(FINALIZE_H5AD.out.versions)

    //
    // Summary report
    //
    SCANPY_GENERATE_REPORT(ch_h5ad, "${projectDir}/modules/local/scanpy/report/templates/scdownstream_analysis_report.qmd", clustering_name, params.report_table_row_limit, params.umap_for_plots, params.umap_color_by ?: '', params.umap_color_by_embeddings ?: '', params.report_max_celltype_labels)


    emit:
    h5ad = ch_h5ad
    versions       = ch_versions // channel: [ path(versions.yml) ]
}
