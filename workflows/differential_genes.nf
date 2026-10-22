/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT MODULES / SUBWORKFLOWS / FUNCTIONS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { DECOUPLER_PSEUDOBULK              } from '../modules/local/decoupler/pseudobulk'
include { FILTER_PSEUDOBULK                 } from '../modules/local/decoupler/filter_pseudobulk'
include { ADATA_SPLITCOL as SPLIT_PER_GROUP } from '../modules/local/adata/splitcol'
include { DIFFERENTIAL_GENES_PER_CONTRAST   } from '../modules/local/pydeseq2/differential_genes'
include { PYDESEQ2_GENERATE_REPORT          } from '../modules/local/pydeseq2/report'
include { PYDESEQ2_EXPORT_MARKERS           } from '../modules/local/pydeseq2/export_markers'
include { softwareVersionsToYAML            } from '../subworkflows/nf-core/utils_nfcore_pipeline'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN MAIN WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow DIFFERENTIAL_GENES {
    take:
    ch_base        // channel: [ val(meta), path(h5ad) ]
    ch_contrasts   // channel: [ val(contrast) ] - map with keys contrast_name, variable, reference_group, target_group, list_of_blocked_variables, exclude_samples_col, exclude_samples_values
    group_col      // string: obs column holding the group label (cluster / cell type)
    sample_col     // string: obs column identifying the biological sample
    min_counts     // integer
    min_cells      // integer
    min_samples    // integer

    main:

    ch_versions = Channel.empty()

    ch_h5ad = ch_base

    // Every column that has to survive the pseudobulk aggregation: the group label plus
    // the variable and the blocking variables of every contrast. DECOUPLER_PSEUDOBULK
    // expects a single comma-separated string because its template splits on ','; the
    // report reuses the same string to decide which columns to break the profiles down by.
    // It is a value channel, so both consumers can read it.
    ch_group_cols = ch_contrasts
        .flatMap { contrast ->
            [contrast.variable] + (contrast.list_of_blocked_variables ? contrast.list_of_blocked_variables.split(',').collect { it.trim() } : [])
        }
        .collect()
        .map { cols -> ([group_col] + cols).findAll { it }.unique().join(',') }

    //
    // MODULE: Pseudobulk per (sample x group label x contrast variables)
    //
    DECOUPLER_PSEUDOBULK(ch_h5ad, sample_col, ch_group_cols)
    ch_versions = ch_versions.mix(DECOUPLER_PSEUDOBULK.out.versions)

    //
    // MODULE: Drop pseudobulk samples with too few counts or cells
    //
    FILTER_PSEUDOBULK(DECOUPLER_PSEUDOBULK.out.h5ad, min_counts, min_cells)
    ch_versions = ch_versions.mix(FILTER_PSEUDOBULK.out.versions)

    //
    // MODULE: Split the pseudobulk object into one object per instance of the group label
    //
    SPLIT_PER_GROUP(FILTER_PSEUDOBULK.out.h5ad, group_col)
    ch_versions = ch_versions.mix(SPLIT_PER_GROUP.out.versions)

    // The split emits a list of h5ads, one per group instance - flatten it and recover the
    // group value from the file name. Each subset is then tested for every contrast, with
    // both facets carried in the meta so the outputs stay distinguishable.
    ch_de_input = SPLIT_PER_GROUP.out.h5ad
        .transpose()
        .map { meta, h5ad -> [meta + [subset: h5ad.simpleName], h5ad] }
        .combine(ch_contrasts)
        .map { meta, h5ad, contrast -> [meta + [contrast: contrast.contrast_name], h5ad, contrast] }

    //
    // MODULE: Differential genes for one contrast within one instance of the group label
    //
    DIFFERENTIAL_GENES_PER_CONTRAST(ch_de_input, min_samples)
    ch_versions = ch_versions.mix(DIFFERENTIAL_GENES_PER_CONTRAST.out.versions)

    // Fan-in: the report summarises every result at once, so the per-(group instance x
    // contrast) channel is collapsed to a single collection of TSVs.
    ch_de_tsvs = DIFFERENTIAL_GENES_PER_CONTRAST.out.tsv
        .map { _meta, tsv -> tsv }
        .collect()

    // The group instance and the contrast name are recorded next to each file name instead
    // of being parsed back out of it - both ids may contain underscores, so the file name
    // alone is ambiguous. Headerless, one row per result file, deterministically sorted.
    ch_de_manifest = DIFFERENTIAL_GENES_PER_CONTRAST.out.tsv
        .map { meta, tsv -> "${tsv.name}\t${meta.subset}\t${meta.contrast}" }
        .collectFile(name: 'de_manifest.tsv', newLine: true, sort: true)

    // The marker export additionally needs the contrast definition (variable, target and
    // reference level) next to each result file. It gets a manifest of its own so that meta
    // and the report's ch_de_manifest stay as they are. The contrast is looked up by name
    // with combine(by: 0) rather than join: every contrast yields one result per group
    // instance, and join would pair only the first of them with the contrast. meta.contrast
    // is taken from ch_contrasts upstream, so every result finds its contrast.
    ch_de_markers_manifest = DIFFERENTIAL_GENES_PER_CONTRAST.out.tsv
        .map { meta, tsv -> [meta.contrast, meta.subset, tsv] }
        .combine(ch_contrasts.map { contrast -> [contrast.contrast_name, contrast] }, by: 0)
        .map { contrast_name, subset, tsv, contrast ->
            "${tsv.name}\t${subset}\t${contrast_name}\t${contrast.variable}\t${contrast.target_group}\t${contrast.reference_group}"
        }
        .collectFile(name: 'de_markers_manifest.tsv', newLine: true, sort: true)

    // The contrast definitions, so the report can state what each contrast compares. Same
    // columns as the contrasts file, one row per contrast, deterministically sorted below
    // the header. Blocking variables and excluded values stay comma-separated as given.
    ch_contrasts_table = ch_contrasts
        .map { contrast ->
            [
                contrast.contrast_name,
                contrast.variable,
                contrast.reference_group,
                contrast.target_group,
                contrast.list_of_blocked_variables,
                contrast.exclude_samples_col,
                contrast.exclude_samples_values,
            ].collect { value -> value ?: '' }.join('\t')
        }
        .collectFile(
            name: 'contrasts.tsv',
            seed: "id\tvariable\treference\ttarget\tblocking\texclude_samples_col\texclude_samples_values",
            newLine: true,
            sort: true,
        )

    // Versions of the modules that produced the results above, collated into a single file
    // so the report can list them. The report's own versions are mixed into ch_versions
    // only afterwards - feeding them back in here would close a cycle.
    ch_module_versions = softwareVersionsToYAML(ch_versions)
        .collectFile(name: 'differential_genes_versions.yml', sort: true, newLine: true)

    //
    // MODULE: Quarto report over all differential expression results
    //
    PYDESEQ2_GENERATE_REPORT(
        FILTER_PSEUDOBULK.out.h5ad,
        ch_de_tsvs,
        ch_de_manifest,
        ch_contrasts_table,
        "${projectDir}/modules/local/pydeseq2/report/templates/scdownstream_differential_genes_report.qmd",
        sample_col,
        ch_group_cols,
        params.report_table_row_limit,
        ch_module_versions,
    )
    ch_versions = ch_versions.mix(PYDESEQ2_GENERATE_REPORT.out.versions)

    //
    // MODULE: Marker-gene JSON with the top differential genes per group label, contrast and direction
    //
    // Its versions join ch_versions only here, after ch_module_versions was collated for the
    // report, so the report does not depend on this module.
    PYDESEQ2_EXPORT_MARKERS(
        ch_de_tsvs,
        ch_de_markers_manifest,
        params.name ?: 'markers',
        group_col,
        params.diffgenes_markers_thr_adj_pvalue,
        params.diffgenes_markers_n_top,
        params.diffgenes_markers_min_logfc,
        params.markers_sep,
    )
    ch_versions = ch_versions.mix(PYDESEQ2_EXPORT_MARKERS.out.versions)

    //
    // Collate and save software versions
    //
    softwareVersionsToYAML(ch_versions)
        .collectFile(
            storeDir: "${params.outdir}/pipeline_info",
            name: 'nf_core_' + 'scdownstream_differential_genes_software_' + 'mqc_' + 'versions.yml',
            sort: true,
            newLine: true,
        )

    emit:
    tsv             = DIFFERENTIAL_GENES_PER_CONTRAST.out.tsv // channel: [ val(meta), path(tsv) ]
    h5ad_pseudobulk = FILTER_PSEUDOBULK.out.h5ad              // channel: [ val(meta), path(h5ad) ]
    html            = PYDESEQ2_GENERATE_REPORT.out.html       // channel: [ val(meta), path(html) ]
    versions        = ch_versions                             // channel: [ path(versions.yml) ]
}
