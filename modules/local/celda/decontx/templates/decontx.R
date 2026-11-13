#!/usr/bin/env Rscript

library(anndataR)
library(celda)

# Read the AnnData object
adata <- read_h5ad("${h5ad}")

# Convert to SingleCellExperiment
sce <- adata\$as_SingleCellExperiment()

# Prepare parameters for decontX
# anndataR 1.2.1 converts X to the assay "X" and each layer to an assay of the same name, so the assay is
# the input layer itself
params <- list()
params\$assayName <- "${input_layer}"

# Handle batch information if available
batch_col <- "${batch_col}"
if (batch_col != "" && length(unique(adata\$obs[[batch_col]])) > 1) {
    params\$batch <- adata\$obs[[batch_col]]
}

# Handle background data if available
raw_path <- "${raw}"
if (file.exists(raw_path)) {
    raw <- read_h5ad(raw_path)
    params\$background <- raw\$as_SingleCellExperiment()
    # The raw counts are in the same layer as in the filtered object, or else in X
    bg_assays <- SummarizedExperiment::assayNames(params\$background)
    params\$bgAssayName <- if (params\$assayName %in% bg_assays) params\$assayName else "X"
    message("decontX background: raw counts from ", ifelse(params\$bgAssayName == "X", "X", paste0("layers['", params\$bgAssayName, "']")))
}

# Run decontX with parameters
corrected <- do.call(decontX, c(list(sce), params))

# The corrected counts as cells x genes. Taken from the assay directly, not through as_AnnData(corrected),
# which in anndataR 1.2.1 would also convert the input assay "X" and decontX's results in metadata()
decontx_counts <- Matrix::t(decontXcounts(corrected))

# Update X or the output layer
if ("${output_layer}" == "X") {
    adata\$X <- decontx_counts
} else {
    adata\$layers[["${output_layer}"]] <- decontx_counts
}

# Save the output
write_h5ad(adata, "${prefix}.h5ad")

# Write version information
writeLines(
    c(
        '"${task.process}":',
        paste('    r:', paste(version\$major, version\$minor, sep = ".")),
        paste('    anndataR:', as.character(packageVersion('anndataR'))),
        paste('    celda:', as.character(packageVersion('celda')))
    ),
'versions.yml')
