#!/usr/bin/env nextflow
/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    nf-core/scdownstream
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    Github : https://github.com/nf-core/scdownstream
    Website: https://nf-co.re/scdownstream
    Slack  : https://nfcore.slack.com/channels/scdownstream
----------------------------------------------------------------------------------------
*/

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT FUNCTIONS / MODULES / SUBWORKFLOWS / WORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { QC_CLUSTER              } from './workflows/qc_clustering'
include { DOWNSTREAM_ANALYSIS     } from './workflows/downstream'
include { DIFFERENTIAL_GENES      } from './workflows/differential_genes'
include { PIPELINE_INITIALISATION } from './subworkflows/local/utils_nfcore_scdownstream_pipeline'
include { PIPELINE_COMPLETION     } from './subworkflows/local/utils_nfcore_scdownstream_pipeline'
include { samplesheetToList       } from 'plugin/nf-schema'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN QC and CLUSTERING WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
workflow qc_clustering {

    main:
    //
    // SUBWORKFLOW: Run initialisation tasks
    //
    PIPELINE_INITIALISATION (
        params.version,
        params.validate_params,
        params.monochrome_logs,
        args,
        params.outdir,
        params.input
    )

    //
    // WORKFLOW: Run main workflow
    //
    QC_CLUSTER (
        PIPELINE_INITIALISATION.out.samplesheet,
        params.base_adata
            ? Channel.value([[id: params.name ? params.name : "qc_clustering" ], file(params.base_adata, checkIfExists: true)])
            : Channel.value([[], []])
    )

    //
    // SUBWORKFLOW: Run completion tasks
    //
    PIPELINE_COMPLETION (
        params.email,
        params.email_on_fail,
        params.plaintext_email,
        params.outdir,
        params.monochrome_logs,
        params.hook_url,
        QC_CLUSTER.out.multiqc_report
    )
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN DOWNSTREAM WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
workflow downstream {

    main:
    //
    // SUBWORKFLOW: Run initialisation tasks
    //
    PIPELINE_INITIALISATION (
        params.version,
        params.validate_params,
        params.monochrome_logs,
        args,
        params.outdir,
        params.input
    )

    //
    // WORKFLOW: Run main workflow
    //
    DOWNSTREAM_ANALYSIS (
        PIPELINE_INITIALISATION.out.samplesheet,
        params.base_adata
            ? Channel.value([[id: params.name ? params.name : "downstream" ], file(params.base_adata, checkIfExists: true)])
            : Channel.value([[], []])
    )

}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN DIFFERENTIAL GENES WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
workflow differential_genes {

    main:
    //
    // SUBWORKFLOW: Run initialisation tasks
    //
    PIPELINE_INITIALISATION (
        params.version,
        params.validate_params,
        params.monochrome_logs,
        args,
        params.outdir,
        params.input
    )

    //
    // Read the contrasts file following the nf-core/differentialabundance definition
    //
    ch_contrasts = Channel.fromList(samplesheetToList(params.diffgenes_contrasts, "${projectDir}/assets/schema_contrasts.json"))
        .map { id, variable, reference, target, blocking, exclude_samples_col, exclude_samples_values ->
            [
                contrast_name: id,
                variable: variable,
                reference_group: reference,
                target_group: target,
                // the contrasts file is tab-separated, so blocking variables are comma-separated,
                // which is exactly what the DE template splits on
                list_of_blocked_variables: blocking ?: '',
                exclude_samples_col: exclude_samples_col ?: '',
                exclude_samples_values: exclude_samples_values ?: ''
            ]
        }

    //
    // WORKFLOW: Run main workflow
    //
    DIFFERENTIAL_GENES (
        params.base_adata
            ? Channel.value([[id: params.name ? params.name : "differential_genes" ], file(params.base_adata, checkIfExists: true)])
            : Channel.value([[], []]),
        ch_contrasts,
        params.diffgenes_group_col,
        params.diffgenes_sample_col,
        params.diffgenes_min_counts,
        params.diffgenes_min_cells,
        params.diffgenes_min_samples
    )

}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN MAIN WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow {

    main:
    // Without -entry, run the first stage. Users should always choose the stage explicitly.
    log.warn "No -entry given: running qc_clustering. Pass -entry qc_clustering, downstream or differential_genes explicitly."
    qc_clustering()
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    THE END
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
