/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT MODULES / SUBWORKFLOWS / FUNCTIONS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/


include { MULTIQC                              } from '../modules/nf-core/multiqc/main'
include { paramsSummaryMap                     } from 'plugin/nf-schema'
include { paramsSummaryMultiqc                 } from '../subworkflows/nf-core/utils_nfcore_pipeline'
include { softwareVersionsToYAML               } from '../subworkflows/nf-core/utils_nfcore_pipeline'
include { methodsDescriptionText               } from '../subworkflows/local/utils_nfcore_scdownstream_pipeline'
include { SCANPY_RANKGENESGROUPS               } from '../modules/local/scanpy/rankgenesgroups'
include { SCANPY_ENRICH                        } from '../modules/local/scanpy/enrich'
include { SCANPY_GENERATE_REPORT               } from '../modules/local/scanpy/report'
include { SCANPY_REPORT_TO_HTML                } from '../modules/local/scanpy/report'
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

    SCANPY_RANKGENESGROUPS(ch_h5ad, clustering_name)
    ch_h5ad = SCANPY_RANKGENESGROUPS.out.h5ad
    SCANPY_ENRICH(ch_h5ad, "rank_genes_groups", params.species, params.enrich_min_in_group_fraction, params.enrich_min_fold_change, params.enrich_max_out_group_fraction)
    ch_h5ad = SCANPY_ENRICH.out.h5ad
    LIANA_RANKAGGREGATE(ch_h5ad, clustering_name, params.species, params.ortholog_hcop_directory)
    ch_h5ad = LIANA_RANKAGGREGATE.out.h5ad

    // create output final output files
    FINALIZE_H5AD(ch_h5ad)
    
    ch_versions = ch_versions.mix(FINALIZE_H5AD.out.versions)

    //
    // Summary report
    //
    SCANPY_GENERATE_REPORT(ch_h5ad, "${projectDir}/modules/local/scanpy/report/templates/scdownstream_analysis_report.ipynb", clustering_name, params.report_table_row_limit)
    SCANPY_REPORT_TO_HTML(SCANPY_GENERATE_REPORT.out.ipynp)


    emit:
    h5ad = ch_h5ad
    versions       = ch_versions // channel: [ path(versions.yml) ]
}
