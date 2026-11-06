//
// Subworkflow with functionality specific to the nf-core/scdownstream pipeline
//

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT FUNCTIONS / MODULES / SUBWORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { UTILS_NFSCHEMA_PLUGIN   } from '../../nf-core/utils_nfschema_plugin'
include { paramsSummaryMap        } from 'plugin/nf-schema'
include { samplesheetToList       } from 'plugin/nf-schema'
include { completionEmail         } from '../../nf-core/utils_nfcore_pipeline'
include { completionSummary       } from '../../nf-core/utils_nfcore_pipeline'
include { imNotification          } from '../../nf-core/utils_nfcore_pipeline'
include { UTILS_NFCORE_PIPELINE   } from '../../nf-core/utils_nfcore_pipeline'
include { UTILS_NEXTFLOW_PIPELINE } from '../../nf-core/utils_nextflow_pipeline'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    SUBWORKFLOW TO INITIALISE PIPELINE
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow PIPELINE_INITIALISATION {
    take:
    version           // boolean: Display version and exit
    validate_params   // boolean: Boolean whether to validate parameters against the schema at runtime
    monochrome_logs   // boolean: Do not use coloured log outputs
    nextflow_cli_args //   array: List of positional nextflow CLI args
    outdir            //  string: The output directory where the results will be saved
    input             //  string: Path to input samplesheet

    main:

    ch_versions = Channel.empty()

    //
    // Print version and exit if required and dump pipeline parameters to JSON file
    //
    UTILS_NEXTFLOW_PIPELINE(
        version,
        true,
        outdir,
        workflow.profile.tokenize(',').intersect(['conda', 'mamba']).size() >= 1,
    )

    //
    // Validate parameters and generate parameter summary to stdout
    //
    UTILS_NFSCHEMA_PLUGIN(
        workflow,
        validate_params,
        null,
    )

    //
    // Check parameter combinations the schema cannot express
    //
    validatePcaParameters()
    if (params.stage == 'qc_clustering') {
        validateInputParameters()
    }
    if (params.stage == 'train_models') {
        validateTrainParameters()
    }

    //
    // Check config provided to the pipeline
    //
    UTILS_NFCORE_PIPELINE(
        nextflow_cli_args
    )

    //
    // Create channel from input file provided through params.input
    //
    ch_samplesheet = params.input
        ? Channel.fromList(samplesheetToList(params.input, "${projectDir}/assets/schema_input.json")).map {
            validateInputSamplesheet(it)
        }
        : Channel.empty()

    emit:
    samplesheet = ch_samplesheet
    versions    = ch_versions
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    SUBWORKFLOW FOR PIPELINE COMPLETION
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow PIPELINE_COMPLETION {
    take:
    email           //  string: email address
    email_on_fail   //  string: email address sent on pipeline failure
    plaintext_email // boolean: Send plain-text email instead of HTML
    outdir          //    path: Path to output directory where results will be published
    monochrome_logs // boolean: Disable ANSI colour codes in log output
    hook_url        //  string: hook URL for notifications
    multiqc_report  //  string: Path to MultiQC report

    main:
    summary_params = paramsSummaryMap(workflow, parameters_schema: "nextflow_schema.json")
    def multiqc_reports = multiqc_report.toList()

    //
    // Completion email and summary
    //
    workflow.onComplete {
        if (email || email_on_fail) {
            completionEmail(
                summary_params,
                email,
                email_on_fail,
                plaintext_email,
                outdir,
                monochrome_logs,
                multiqc_reports.getVal(),
            )
        }

        completionSummary(monochrome_logs)
        if (hook_url) {
            imNotification(summary_params, hook_url)
        }
    }

    workflow.onError {
        log.error("Pipeline failed. Please refer to troubleshooting docs: https://nf-co.re/docs/usage/troubleshooting")
    }
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    FUNCTIONS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
//
// Check the integration parameter combinations of --stage qc_clustering
//
def validateInputParameters() {
    def integration_methods = params.integration_methods
        ? params.integration_methods.split(',').collect { it -> it.trim().toLowerCase() }
        : []

    if (params.scvi_model && !('scvi' in integration_methods)) {
        throw new Exception("--scvi_model requires scvi in --integration_methods.")
    }

    if (params.scanvi_model && !('scanvi' in integration_methods)) {
        throw new Exception("--scanvi_model requires scanvi in --integration_methods.")
    }

    // INTEGRATE builds its output object from scVI, so scANVI alone would produce nothing
    if ('scanvi' in integration_methods && !('scvi' in integration_methods)) {
        throw new Exception("scanvi in --integration_methods requires scvi as well.")
    }

    if (params.base_adata && 'scanvi' in integration_methods) {
        throw new Exception("Reference mapping onto --base_adata with scANVI is not supported yet.")
    }

    // Extending --base_adata with the --input samples. Without --input, the base is processed
    // on its own (tests/main_pipeline_sub.nf.test) and these extension checks do not apply.
    if (params.input && params.base_adata && (integration_methods - ['scvi', 'scanvi', 'scimilarity']).size() > 0) {
        throw new Exception("Only scvi, scanvi and scimilarity integration methods are supported if base_adata is provided")
    }

    if (params.input && params.base_adata && 'scvi' in integration_methods && !params.scvi_model) {
        throw new Exception("If base_adata is provided and scvi is used as integration method, scvi_model must be provided.")
    }

    if (params.input && params.base_adata && 'scimilarity' in integration_methods && !params.scimilarity_model) {
        throw new Exception("If base_adata is provided and scimilarity is used as integration method, scimilarity_model must be provided.")
    }
}

//
// Read a boolean parameter. With the Nextflow 26 strict parser and no typed params block,
// parameters given on the command line stay strings: `--train_evaluate false` gives "false",
// which Groovy treats as true. A -params-file (YAML/JSON) and nextflow.config keep the type.
// Booleans are returned as they are, null is false, and strings are true only for
// "true", "yes" or "1" (case-insensitive).
//
def paramAsBoolean(value) {
    if (value instanceof Boolean) {
        return value
    }
    if (value == null) {
        return false
    }
    return value.toString().trim().toLowerCase() in ['true', 'yes', '1']
}

//
// Check the parameters of --stage train_models. It reads --train_adata, not --input.
//
def validateTrainParameters() {
    if (!params.train_adata) {
        throw new Exception("--stage train_models requires --train_adata, the annotated reference h5ad.")
    }

    if (!params.train_label_col) {
        throw new Exception("--stage train_models requires --train_label_col, the obs column with the cell type labels.")
    }

    def supported_methods = ['celltypist', 'scanvi']
    def train_methods = params.train_methods
        ? params.train_methods.toString().split(',').collect { it -> it.trim().toLowerCase() }.findAll { it -> it }
        : []
    if (!train_methods) {
        throw new Exception("--train_methods must name at least one method (${supported_methods.join(', ')}).")
    }
    def unsupported_methods = train_methods - supported_methods
    if (unsupported_methods) {
        throw new Exception("--train_methods: unsupported method(s) ${unsupported_methods.join(', ')}. Supported: ${supported_methods.join(', ')}.")
    }

    // Numbers given on the command line arrive as strings (see paramAsBoolean), so convert
    // them before comparing
    def test_fraction = params.train_test_fraction as BigDecimal
    if (!(test_fraction > 0 && test_fraction < 1)) {
        throw new Exception("--train_test_fraction (${params.train_test_fraction}) must be larger than 0 and smaller than 1.")
    }

    def min_cells_per_label = params.train_min_cells_per_label as Integer
    if (!(min_cells_per_label >= 1)) {
        throw new Exception("--train_min_cells_per_label (${params.train_min_cells_per_label}) must be at least 1.")
    }

    if (params.train_max_cells_per_label != null && !((params.train_max_cells_per_label as Integer) >= 1)) {
        throw new Exception("--train_max_cells_per_label (${params.train_max_cells_per_label}) must be at least 1.")
    }

    if (paramAsBoolean(params.train_celltypist_mini_batch) && !paramAsBoolean(params.train_celltypist_use_sgd)) {
        throw new Exception("--train_celltypist_mini_batch requires --train_celltypist_use_sgd.")
    }

    // Without the evaluation there is no training-split model to keep, so the final models are
    // always trained on all cells
    if (!paramAsBoolean(params.train_evaluate) && !paramAsBoolean(params.train_final_on_all)) {
        log.warn("--train_final_on_all false is ignored with --train_evaluate false: the final models are trained on all cells.")
    }
}

//
// Check that the PCA neighbour graph only asks for PCs that SCANPY_PCA computes
//
def validatePcaParameters() {
    if (params.neighbors_n_pcs && params.pca_n_comps && params.neighbors_n_pcs > params.pca_n_comps) {
        throw new Exception("--neighbors_n_pcs (${params.neighbors_n_pcs}) cannot be larger than --pca_n_comps (${params.pca_n_comps}): the PCA graph can only use PCs that were computed.")
    }
}

//
// Validate channels from input samplesheet
//
def validateInputSamplesheet(input) {
    def (meta, filtered, unfiltered) = input
    if (!filtered && !unfiltered) {
        throw new Exception("Both filtered and unfiltered files are missing for sample ${meta.id}")
    }

    return input
}

//
// Return the per-sample value meta[key] from the samplesheet, or `global` when the
// samplesheet cell was left empty (null, '' or []). Explicit 0 and false are kept.
//
def sampleOrGlobal(meta, key, global) {
    def value = meta[key]
    if (value == null || value == '' || value == []) {
        return global
    }
    return value
}

//
// Generate methods description for MultiQC
//
def toolCitationText() {
    // TODO nf-core: Optionally add in-text citation tools to this list.
    // Can use ternary operators to dynamically construct based conditions, e.g. params["run_xyz"] ? "Tool (Foo et al. 2023)" : "",
    // Uncomment function in methodsDescriptionText to render in MultiQC report
    def citation_text = [
        "Tools used in the workflow included:",
        "MultiQC (Ewels et al. 2016)",
        ".",
    ].join(' ').trim()

    return citation_text
}

def toolBibliographyText() {
    // TODO nf-core: Optionally add bibliographic entries to this list.
    // Can use ternary operators to dynamically construct based conditions, e.g. params["run_xyz"] ? "<li>Author (2023) Pub name, Journal, DOI</li>" : "",
    // Uncomment function in methodsDescriptionText to render in MultiQC report
    def reference_text = [
        "<li>Ewels, P., Magnusson, M., Lundin, S., & Käller, M. (2016). MultiQC: summarize analysis results for multiple tools and samples in a single report. Bioinformatics , 32(19), 3047–3048. doi: /10.1093/bioinformatics/btw354</li>"
    ].join(' ').trim()

    return reference_text
}

def methodsDescriptionText(mqc_methods_yaml) {
    // Convert  to a named map so can be used as with familiar NXF ${workflow} variable syntax in the MultiQC YML file
    def meta = [:]
    meta.workflow = workflow.toMap()
    meta["manifest_map"] = workflow.manifest.toMap()

    // Pipeline DOI
    if (meta.manifest_map.doi) {
        // Using a loop to handle multiple DOIs
        // Removing `https://doi.org/` to handle pipelines using DOIs vs DOI resolvers
        // Removing ` ` since the manifest.doi is a string and not a proper list
        def temp_doi_ref = ""
        def manifest_doi = meta.manifest_map.doi.tokenize(",")
        manifest_doi.each { doi_ref ->
            temp_doi_ref += "(doi: <a href=\'https://doi.org/${doi_ref.replace("https://doi.org/", "").replace(" ", "")}\'>${doi_ref.replace("https://doi.org/", "").replace(" ", "")}</a>), "
        }
        meta["doi_text"] = temp_doi_ref.substring(0, temp_doi_ref.length() - 2)
    }
    else {
        meta["doi_text"] = ""
    }
    meta["nodoi_text"] = meta.manifest_map.doi ? "" : "<li>If available, make sure to update the text to include the Zenodo DOI of version of the pipeline used. </li>"

    // Tool references
    meta["tool_citations"] = ""
    meta["tool_bibliography"] = ""

    // TODO nf-core: Only uncomment below if logic in toolCitationText/toolBibliographyText has been filled!
    // meta["tool_citations"] = toolCitationText().replaceAll(", \\.", ".").replaceAll("\\. \\.", ".").replaceAll(", \\.", ".")
    // meta["tool_bibliography"] = toolBibliographyText()


    def methods_text = mqc_methods_yaml.text

    def engine = new groovy.text.SimpleTemplateEngine()
    def description_html = engine.createTemplate(methods_text).make(meta)

    return description_html.toString()
}
