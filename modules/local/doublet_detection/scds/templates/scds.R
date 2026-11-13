#!/usr/bin/env Rscript

library(scds)
library(SingleCellExperiment)
library(anndataR)

Sys.setenv(TMP = ".")

set.seed(0)

adata <- read_h5ad("${h5ad}")
# scds reads the assay "counts". anndataR 1.2.1 converts X to the assay "X" and each layer to an assay of the
# same name, so without a counts layer X (which holds the raw counts before normalisation) becomes "counts".
if ("counts" %in% adata\$layers_keys()) {
    message("scds input: raw counts from layers['counts']")
    sce <- adata\$as_SingleCellExperiment()
} else {
    message("scds input: raw counts from X")
    sce <- adata\$as_SingleCellExperiment(x_mapping = "counts")
}

## Annotate doublet using binary classification based doublet scoring:
sce <- bcds(sce, retRes = TRUE, estNdbl=TRUE)

## Annotate doublet using co-expression based doublet scoring:
try({
    sce <- cxds(sce, retRes = TRUE, estNdbl=TRUE)
})

### If cxds worked, run hybrid, otherwise use bcds annotations
if ("cxds_score" %in% colnames(colData(sce))) {
    ## Combine both annotations into a hybrid annotation
    sce <- cxds_bcds_hybrid(sce, estNdbl=TRUE)

    predictions <- colData(sce)[, 'hybrid_call', drop=FALSE]
} else {
    predictions <- colData(sce)[, 'bcds_call', drop=FALSE]
}

# Add the scds scores and calls to the input object. Not as_AnnData(sce): anndataR 1.2.1 copies metadata(sce)
# to uns, and with retRes = TRUE that holds the xgboost model and xgb.cv result, which write_h5ad() cannot write.
scds_cols <- grep("^(bcds|cxds|hybrid)_(score|call)\$", colnames(colData(sce)), value = TRUE)
obs <- adata\$obs
obs[scds_cols] <- as.data.frame(colData(sce)[, scds_cols, drop = FALSE])
adata\$obs <- obs
write_h5ad(adata, "${prefix}.h5ad")

colnames(predictions) <- "${prefix}"
write.csv(predictions, "${prefix}.csv")

################################################
################################################
## VERSIONS FILE                              ##
################################################
################################################

r.version <- strsplit(version[['version.string']], ' ')[[1]][3]
scds.version <- as.character(packageVersion('scds'))

writeLines(
    c(
        '"${task.process}":',
        paste('    R:', r.version),
        paste('    scds:', scds.version)
    ),
'versions.yml')

################################################
################################################
################################################
################################################
