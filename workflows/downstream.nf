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
    ch_multiqc_files = Channel.empty()

    //
    // Perform downstream analysis
    //
    clustering_name = params.selected_clustering == null ? 'leiden' :  params.selected_clustering

    SCANPY_RANKGENESGROUPS(ch_base, clustering_name)
    ch_h5ad = SCANPY_RANKGENESGROUPS.out.h5ad
    SCANPY_ENRICH(ch_h5ad, "rank_genes_groups", params.species, params.enrich_min_in_group_fraction, params.enrich_min_fold_change, params.enrich_max_out_group_fraction)
    ch_h5ad = SCANPY_ENRICH.out.h5ad
    LIANA_RANKAGGREGATE(ch_h5ad, params.species, clustering_name)
    ch_h5ad = LIANA_RANKAGGREGATE.out.h5ad

    //
    // Summary report
    //
    SCANPY_GENERATE_REPORT(ch_h5ad, "${projectDir}/modules/local/scanpy/report/templates/scdownstream_analysis_report.ipynb", clustering_name)
    SCANPY_REPORT_TO_HTML(SCANPY_GENERATE_REPORT.out.ipynp)


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
