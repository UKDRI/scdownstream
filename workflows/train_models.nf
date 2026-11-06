/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT MODULES / SUBWORKFLOWS / FUNCTIONS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { ADATA_TRAINSPLIT                                } from '../modules/local/adata/trainsplit'
include { CELLTYPIST_TRAIN as CELLTYPIST_TRAIN_EVAL       } from '../modules/local/celltypes/celltypist_train'
include { CELLTYPIST_TRAIN as CELLTYPIST_TRAIN_FINAL      } from '../modules/local/celltypes/celltypist_train'
include { CELLTYPES_CELLTYPIST as CELLTYPIST_PREDICT_TEST } from '../modules/local/celltypes/celltypist'
include { SCANPY_HVGS as HVGS_EVAL                        } from '../modules/local/scanpy/hvgs'
include { SCANPY_HVGS as HVGS_FINAL                       } from '../modules/local/scanpy/hvgs'
include { SCVITOOLS_SCVI as SCVI_EVAL                     } from '../modules/local/scvitools/scvi'
include { SCVITOOLS_SCVI as SCVI_FINAL                    } from '../modules/local/scvitools/scvi'
include { SCVITOOLS_SCANVI as SCANVI_EVAL                 } from '../modules/local/scvitools/scanvi'
include { SCVITOOLS_SCANVI as SCANVI_QUERY_TEST           } from '../modules/local/scvitools/scanvi'
include { SCVITOOLS_SCANVI as SCANVI_FINAL                } from '../modules/local/scvitools/scanvi'
include { MODEL_EVALUATE as EVALUATE_CELLTYPIST           } from '../modules/local/celltypes/evaluate'
include { MODEL_EVALUATE as EVALUATE_SCANVI               } from '../modules/local/celltypes/evaluate'
include { SCANPY_GENERATE_REPORT_TRAINING                 } from '../modules/local/scanpy/report'
include { softwareVersionsToYAML                          } from '../subworkflows/nf-core/utils_nfcore_pipeline'
include { paramAsBoolean                                  } from '../subworkflows/local/utils_nfcore_scdownstream_pipeline'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN MAIN WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow TRAIN_MODELS {
    take:
    ch_train   // channel: [ val(meta), path(h5ad) ] - the annotated reference, meta.id is the model name
    model_name // string: name of the trained models

    main:

    ch_versions = Channel.empty()
    ch_evaluation = Channel.empty()
    ch_celltypist_model = Channel.empty()
    ch_scanvi_model = Channel.empty()
    ch_scvi_model = Channel.empty()
    ch_reference_h5ad = Channel.empty()

    methods = params.train_methods.split(',').collect { method -> method.trim().toLowerCase() }
    // Booleans given on the command line arrive as strings ("false" is truthy), see paramAsBoolean
    evaluate = paramAsBoolean(params.train_evaluate)
    // Without an evaluation there is no training-split model to keep, so the final models are
    // always fitted on all cells (validateTrainParameters() warns about the ignored setting)
    final_on_all = !evaluate || paramAsBoolean(params.train_final_on_all)

    // The models trained on the training split are named <name>_eval, unless they are kept as
    // the final models (--train_final_on_all false): then they carry the model name itself
    eval_name = final_on_all ? model_name + '_eval' : model_name

    // Every derived meta keeps the model name, so the module config can name the published
    // models after it independently of the task ids below
    ch_train = ch_train.map { meta, h5ad -> [meta + [name: model_name], h5ad] }

    //
    // MODULE: Prepare the reference (counts, labels, batch, symbols), drop small labels, split
    //
    ADATA_TRAINSPLIT(ch_train)
    ch_versions = ch_versions.mix(ADATA_TRAINSPLIT.out.versions)

    //
    // CellTypist
    //
    if (methods.contains('celltypist')) {
        if (evaluate) {
            //
            // MODULE: CellTypist model trained on the training split
            //
            CELLTYPIST_TRAIN_EVAL(
                ADATA_TRAINSPLIT.out.train.map { meta, h5ad -> [meta + [id: eval_name], h5ad] }
            )
            ch_versions = ch_versions.mix(CELLTYPIST_TRAIN_EVAL.out.versions)

            //
            // MODULE: Predict the test cells with the training-split model (a local .pkl, no built-in names)
            //
            CELLTYPIST_PREDICT_TEST(
                ADATA_TRAINSPLIT.out.test.map { meta, h5ad -> [meta + [id: eval_name + '_celltypist_test'], h5ad, 'index'] },
                [],
                CELLTYPIST_TRAIN_EVAL.out.model.map { _meta, pkl -> pkl },
            )
            ch_versions = ch_versions.mix(CELLTYPIST_PREDICT_TEST.out.versions)

            //
            // MODULE: Compare the predictions with the hidden true labels
            //
            // CELLTYPES_CELLTYPIST names the prediction column after the model file stem, so the
            // column is taken from the file itself rather than from the expected name.
            EVALUATE_CELLTYPIST(
                CELLTYPIST_PREDICT_TEST.out.h5ad.map { meta, h5ad -> [meta + [id: eval_name + '_celltypist', method: 'celltypist'], h5ad] },
                CELLTYPIST_TRAIN_EVAL.out.model.map { _meta, pkl -> "celltypist:${pkl.baseName}" },
                'label_true',
            )
            ch_versions = ch_versions.mix(EVALUATE_CELLTYPIST.out.versions)
            ch_evaluation = ch_evaluation.mix(
                EVALUATE_CELLTYPIST.out.metrics,
                EVALUATE_CELLTYPIST.out.per_label,
                EVALUATE_CELLTYPIST.out.confusion,
            )
            ch_celltypist_model = CELLTYPIST_TRAIN_EVAL.out.model
        }

        if (final_on_all) {
            //
            // MODULE: Final CellTypist model, trained on all cells
            //
            CELLTYPIST_TRAIN_FINAL(
                ADATA_TRAINSPLIT.out.all.map { meta, h5ad -> [meta + [id: model_name], h5ad] }
            )
            ch_versions = ch_versions.mix(CELLTYPIST_TRAIN_FINAL.out.versions)
            ch_celltypist_model = CELLTYPIST_TRAIN_FINAL.out.model
        }
    }

    //
    // scVI -> scANVI
    //
    if (methods.contains('scanvi')) {
        if (evaluate) {
            //
            // MODULE: Reference genes of the training split (Pearson residuals on counts, per batch; see conf/modules.config)
            //
            HVGS_EVAL(
                ADATA_TRAINSPLIT.out.train.map { meta, h5ad -> [meta + [id: eval_name], h5ad] },
                params.train_hvgs,
                true,
            )
            ch_versions = ch_versions.mix(HVGS_EVAL.out.versions)

            //
            // MODULE: scVI reference model on the training split
            //
            SCVI_EVAL(
                HVGS_EVAL.out.h5ad.map { meta, h5ad -> [meta + [id: eval_name + '_scvi'], h5ad] },
                [[], []],
                'scratch',
                'batch',
                '',
                '',
            )
            ch_versions = ch_versions.mix(SCVI_EVAL.out.versions)

            //
            // MODULE: scANVI reference model, initialised from the scVI model
            //
            SCANVI_EVAL(
                HVGS_EVAL.out.h5ad.map { meta, h5ad -> [meta + [id: eval_name + '_scanvi'], h5ad] },
                SCVI_EVAL.out.model,
                'from_scvi',
                ['label', 'unknown'],
                'batch',
                '',
                '',
                eval_name,
            )
            ch_versions = ch_versions.mix(SCANVI_EVAL.out.versions)

            //
            // MODULE: Map the test cells (all genes, labels hidden) onto the scANVI reference
            //
            SCANVI_QUERY_TEST(
                ADATA_TRAINSPLIT.out.test.map { meta, h5ad -> [meta + [id: eval_name + '_scanvi_test'], h5ad] },
                SCANVI_EVAL.out.model,
                'query',
                ['label', 'unknown'],
                'batch',
                '',
                '',
                eval_name,
            )
            ch_versions = ch_versions.mix(SCANVI_QUERY_TEST.out.versions)

            //
            // MODULE: Compare the transferred labels with the hidden true labels
            //
            EVALUATE_SCANVI(
                SCANVI_QUERY_TEST.out.h5ad.map { meta, h5ad -> [meta + [id: eval_name + '_scanvi', method: 'scanvi'], h5ad] },
                'scanvi:' + eval_name,
                'label_true',
            )
            ch_versions = ch_versions.mix(EVALUATE_SCANVI.out.versions)
            ch_evaluation = ch_evaluation.mix(
                EVALUATE_SCANVI.out.metrics,
                EVALUATE_SCANVI.out.per_label,
                EVALUATE_SCANVI.out.confusion,
            )
            ch_scvi_model = SCVI_EVAL.out.model
            ch_scanvi_model = SCANVI_EVAL.out.model
            ch_reference_h5ad = SCANVI_EVAL.out.h5ad
        }

        if (final_on_all) {
            //
            // MODULE: Reference genes of all cells
            //
            HVGS_FINAL(
                ADATA_TRAINSPLIT.out.all.map { meta, h5ad -> [meta + [id: model_name], h5ad] },
                params.train_hvgs,
                true,
            )
            ch_versions = ch_versions.mix(HVGS_FINAL.out.versions)

            //
            // MODULE: Final scVI model, trained on all cells
            //
            SCVI_FINAL(
                HVGS_FINAL.out.h5ad.map { meta, h5ad -> [meta + [id: model_name + '_scvi'], h5ad] },
                [[], []],
                'scratch',
                'batch',
                '',
                '',
            )
            ch_versions = ch_versions.mix(SCVI_FINAL.out.versions)

            //
            // MODULE: Final scANVI model, initialised from the final scVI model
            //
            SCANVI_FINAL(
                HVGS_FINAL.out.h5ad.map { meta, h5ad -> [meta + [id: model_name + '_scanvi'], h5ad] },
                SCVI_FINAL.out.model,
                'from_scvi',
                ['label', 'unknown'],
                'batch',
                '',
                '',
                model_name,
            )
            ch_versions = ch_versions.mix(SCANVI_FINAL.out.versions)
            ch_scvi_model = SCVI_FINAL.out.model
            ch_scanvi_model = SCANVI_FINAL.out.model
            ch_reference_h5ad = SCANVI_FINAL.out.h5ad
        }
    }

    // Versions of the modules that produced the results above, collated into a single file
    // so the report can list them. The report's own versions are mixed into ch_versions
    // only afterwards - feeding them back in here would close a cycle.
    ch_module_versions = softwareVersionsToYAML(ch_versions)
        .collectFile(name: 'train_models_versions.yml', sort: true, newLine: true)

    //
    // MODULE: Quarto report on the label counts and the evaluation of every method
    //
    // Fan-in: all evaluation TSVs of all methods go to the single report. collect() of an
    // empty channel emits nothing, which would stop the report, hence the ifEmpty([]).
    SCANPY_GENERATE_REPORT_TRAINING(
        ADATA_TRAINSPLIT.out.label_counts,
        evaluate ? ch_evaluation.map { _meta, tsv -> tsv }.collect().ifEmpty([]) : [],
        "${projectDir}/modules/local/scanpy/report/templates/scdownstream_training_report.qmd",
        ch_module_versions,
    )
    ch_versions = ch_versions.mix(SCANPY_GENERATE_REPORT_TRAINING.out.versions)

    //
    // Collate and save software versions
    //
    softwareVersionsToYAML(ch_versions)
        .collectFile(
            storeDir: "${params.outdir}/pipeline_info",
            name: 'nf_core_' + 'scdownstream_train_models_software_' + 'mqc_' + 'versions.yml',
            sort: true,
            newLine: true,
        )

    emit:
    celltypist_model = ch_celltypist_model                         // channel: [ val(meta), path(pkl) ]
    scvi_model       = ch_scvi_model                               // channel: [ val(meta), path(model.pt) ]
    scanvi_model     = ch_scanvi_model                             // channel: [ val(meta), path(model.pt) ]
    reference_h5ad   = ch_reference_h5ad                           // channel: [ val(meta), path(h5ad) ]
    evaluation       = ch_evaluation                               // channel: [ val(meta), path(tsv) ]
    html             = SCANPY_GENERATE_REPORT_TRAINING.out.html    // channel: [ val(meta), path(html) ]
    versions         = ch_versions                                 // channel: [ path(versions.yml) ]
}
