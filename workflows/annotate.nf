/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT MODULES / SUBWORKFLOWS / FUNCTIONS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { ADATA_PREPANNOTATE                          } from '../modules/local/adata/prepannotate'
include { CELLTYPES_CELLTYPIST as CELLTYPIST_ANNOTATE } from '../modules/local/celltypes/celltypist'
include { SCVITOOLS_SCANVI as SCANVI_ANNOTATE         } from '../modules/local/scvitools/scanvi'
include { ADATA_ADD_OBS_OBSM as ADD_PREDICTIONS       } from '../modules/local/adata/add_obs_obsm'
include { SCANPY_GENERATE_REPORT_ANNOTATION           } from '../modules/local/scanpy/report'
include { FINALIZE_H5AD_SIMPLE as FINALIZE_ANNOTATED  } from '../subworkflows/local/finalize'
include { softwareVersionsToYAML                      } from '../subworkflows/nf-core/utils_nfcore_pipeline'
include { splitCelltypistModels                       } from '../subworkflows/local/utils'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN MAIN WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow ANNOTATE {
    take:
    ch_base // channel: [ val(meta), path(h5ad) ] - the processed object (raw counts in layers['counts'], or in X), meta.id is the output name

    main:

    ch_versions = Channel.empty()
    ch_predictions = Channel.empty()

    // A method runs when its model is given (validateAnnotateParameters() requires at least one)
    celltypist = splitCelltypistModels(params.celltypist_model)
    run_celltypist = celltypist.names || celltypist.files
    // The scANVI prediction columns are named after the folder holding model.pt, as in INTEGRATE
    scanvi_model_name = params.scanvi_model ? file(params.scanvi_model).parent.name : null

    // This run's prediction columns: celltypist:<model> (the built-in name, or the file name
    // without .pkl - splitCelltypistModels() treats every entry ending in .pkl or containing
    // a / as a file) and scanvi:<model>
    prediction_cols = celltypist.names.collect { name -> "celltypist:${name}" } +
        celltypist.files.collect { model_file -> "celltypist:${model_file.baseName}" } +
        (scanvi_model_name ? ["scanvi:${scanvi_model_name}"] : [])

    // Optional cell subset (validateAnnotateParameters() requires both or neither), described
    // for the report as "<column>: <value>, <value>"
    subset_values = params.annotate_subset_values != null
        ? params.annotate_subset_values.toString().split(',').collect { value -> value.trim() }.findAll { value -> value }
        : []
    subset_description = params.annotate_subset_col
        ? "${params.annotate_subset_col}: ${subset_values.join(', ')}"
        : ''

    //
    // MODULE: Working copy for the tools (raw counts in layers['counts'], X as in the input,
    // gene symbols as var_names, obs['batch'], only the cells of the subset if one is given,
    // without earlier values of this run's prediction columns)
    //
    ADATA_PREPANNOTATE(
        ch_base,
        params.annotate_batch_col,
        params.annotate_symbol_col,
        params.annotate_subset_col ?: '',
        params.annotate_subset_values ?: '',
        prediction_cols.join(','),
    )
    ch_versions = ch_versions.mix(ADATA_PREPANNOTATE.out.versions)

    //
    // MODULE: CellTypist predictions with the built-in and local models
    //
    if (run_celltypist) {
        // The working copy has the gene symbols as var_names
        CELLTYPIST_ANNOTATE(
            ADATA_PREPANNOTATE.out.h5ad.map { meta, h5ad -> [meta, h5ad, 'index'] },
            celltypist.names,
            celltypist.files,
        )
        ch_versions = ch_versions.mix(CELLTYPIST_ANNOTATE.out.versions)
        ch_predictions = ch_predictions.mix(CELLTYPIST_ANNOTATE.out.h5ad.map { _meta, h5ad -> h5ad })
    }

    //
    // MODULE: Map the cells onto the scANVI reference model and transfer its labels
    //
    if (params.scanvi_model) {
        // Covariates come from the model's registry in query mode, hence the empty arguments
        SCANVI_ANNOTATE(
            ADATA_PREPANNOTATE.out.h5ad.map { meta, h5ad -> [meta + [model_name: scanvi_model_name], h5ad] },
            Channel.value([[id: 'scanvi_model'], file(params.scanvi_model, checkIfExists: true)]),
            'query',
            ['label', 'unknown'],
            'batch',
            '',
            '',
            scanvi_model_name,
        )
        ch_versions = ch_versions.mix(SCANVI_ANNOTATE.out.versions)
        ch_predictions = ch_predictions.mix(SCANVI_ANNOTATE.out.h5ad.map { _meta, h5ad -> h5ad })
    }

    //
    // MODULE: Copy the predictions into the unchanged base object
    //
    // Fan-in: the one or two tool outputs become a single list of files next to the base
    // object. The extra [ ] keeps combine() from flattening the list into separate elements.
    // meta.prediction_cols and meta.scanvi_model_name select the obs columns and the obsm key
    // to copy (ext.obs_prefixes / ext.obsm_keys in conf/modules.config).
    ADD_PREDICTIONS(
        ch_base
            .map { meta, h5ad -> [meta + [scanvi_model_name: scanvi_model_name, prediction_cols: prediction_cols.join(',')], h5ad] }
            .combine(ch_predictions.collect().map { files -> [files] })
    )
    ch_versions = ch_versions.mix(ADD_PREDICTIONS.out.versions)

    //
    // SUBWORKFLOW: Publish the annotated object as h5ad
    //
    FINALIZE_ANNOTATED(ADD_PREDICTIONS.out.h5ad)
    ch_versions = ch_versions.mix(FINALIZE_ANNOTATED.out.versions)

    // Versions of the modules that produced the results above, collated into a single file
    // so the report can list them. The report's own versions are mixed into ch_versions
    // only afterwards - feeding them back in here would close a cycle.
    ch_module_versions = softwareVersionsToYAML(ch_versions)
        .collectFile(name: 'annotate_versions.yml', sort: true, newLine: true)

    //
    // MODULE: Quarto report on the new predictions
    //
    SCANPY_GENERATE_REPORT_ANNOTATION(
        ADD_PREDICTIONS.out.h5ad,
        "${projectDir}/modules/local/scanpy/report/templates/scdownstream_annotation_report.qmd",
        prediction_cols.join(','),
        subset_description,
        params.umap_for_plots,
        params.report_max_celltype_labels,
        ch_module_versions,
    )
    ch_versions = ch_versions.mix(SCANPY_GENERATE_REPORT_ANNOTATION.out.versions)

    //
    // Collate and save software versions
    //
    softwareVersionsToYAML(ch_versions)
        .collectFile(
            storeDir: "${params.outdir}/pipeline_info",
            name: 'nf_core_' + 'scdownstream_annotate_software_' + 'mqc_' + 'versions.yml',
            sort: true,
            newLine: true,
        )

    emit:
    h5ad     = ADD_PREDICTIONS.out.h5ad                       // channel: [ val(meta), path(h5ad) ]
    html     = SCANPY_GENERATE_REPORT_ANNOTATION.out.html     // channel: [ val(meta), path(html) ]
    versions = ch_versions                                    // channel: [ path(versions.yml) ]
}
