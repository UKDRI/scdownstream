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
include { TRAIN_MODELS            } from './workflows/train_models'
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
            : Channel.value([[], []]),
        params.metadata
            ? Channel.value(file(params.metadata, checkIfExists: true))
            : Channel.value([])
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
    RUN MODEL TRAINING WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
workflow train_models {

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

    // Name of the trained models and of their prediction columns
    train_name = params.train_name ?: params.name ?: 'reference'

    //
    // WORKFLOW: Run main workflow
    //
    TRAIN_MODELS (
        Channel.value([[id: train_name], file(params.train_adata, checkIfExists: true)]),
        train_name
    )

    //
    // SUBWORKFLOW: Run completion tasks
    //
    // There is no MultiQC report in this stage
    PIPELINE_COMPLETION (
        params.email,
        params.email_on_fail,
        params.plaintext_email,
        params.outdir,
        params.monochrome_logs,
        params.hook_url,
        Channel.empty()
    )
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN MAIN WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow {

    main:
    // Dispatch on --stage. The strict syntax parser does not support -entry, so the
    // four named workflows above are selected here instead.
    if (params.help || params.help_full) {
        // nf-schema already printed the help at pipeline start; run no stage, even with --stage
    }
    else if (params.version) {
        // Prints the version and exits
        PIPELINE_INITIALISATION (
            params.version,
            params.validate_params,
            params.monochrome_logs,
            args,
            params.outdir,
            params.input
        )
    }
    else if (params.stage == 'qc_clustering') {
        qc_clustering()
    }
    else if (params.stage == 'downstream') {
        downstream()
    }
    else if (params.stage == 'differential_genes') {
        differential_genes()
    }
    else if (params.stage == 'train_models') {
        train_models()
    }
    else {
        error "Invalid or missing --stage '${params.stage}'. Valid values: qc_clustering, downstream, differential_genes, train_models."
    }
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    THE END
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
