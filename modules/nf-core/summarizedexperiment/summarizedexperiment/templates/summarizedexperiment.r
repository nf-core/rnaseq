#!/usr/bin/env Rscript

# Written by Lorena Pantano and revised for flexibility in handling assays
# Released under the MIT license.

library(SummarizedExperiment)

#' Flexibly read CSV or TSV files
#'
#' @param file Input file
#' @param header Passed to read.delim()
#' @param row.names Passed to read.delim()
#'
#' @return output Data frame

read_delim_flexible <- function(file, header = TRUE, row.names = NULL, check.names = FALSE, stringsAsFactors = FALSE){

    ext <- tolower(tail(strsplit(basename(file), split = "\\\\.")[[1]], 1))

    if (ext == "tsv" || ext == "txt") {
        separator <- "\\t"
    } else if (ext == "csv") {
        separator <- ","
    } else {
        stop(paste("Unknown separator for", ext))
    }

    read.delim(
        file,
        sep = separator,
        header = header,
        row.names = row.names,
        check.names = check.names,
        stringsAsFactors = stringsAsFactors
    )
}

#' Parse out options from a string without recourse to optparse
#'
#' @param x Long-form argument list like --opt1 val1 --opt2 val2
#'
#' @return named list of options and values similar to optparse

parse_args <- function(x){
    args_list <- unlist(strsplit(x, ' ?--')[[1]])[-1]
    args_vals <- lapply(args_list, function(x) scan(text=x, what='character', quiet = TRUE))

    # Ensure the option vectors are length 2 (key/ value) to catch empty ones
    args_vals <- lapply(args_vals, function(z){ length(z) <- 2; z})

    parsed_args <- structure(lapply(args_vals, function(x) x[2]), names = lapply(args_vals, function(x) x[1]))
    parsed_args[! is.na(parsed_args)]
}

#' Find First Column Containing All Specified Entries
#'
#' This function searches through each column of a given data frame to find the
#' first column that contains all of the specified entries in a vector. If such
#' a column is found, the name of the column is returned. If no column matches,
#' the function throws an error.
#'
#' @param namesVector A character vector containing the names to be matched.
#' @param df A data frame within which to search for the column containing all
#'   names specified in `namesVector`.
#'
#' @return The name of the first column in `df` that contains all entries from
#'   `namesVector`. If no such column exists, the function will throw an error.

findColumnWithAllEntries <- function(namesVector, df) {
    for (colName in names(df)) {
        if (all(namesVector %in% df[[colName]])) {
            return(colName)
        }
    }
    cat(capture.output(print(df)), sep="\n", file=stderr())
    stop(paste("No column contains all vector entries ", paste(namesVector, collapse = ', ')))
}

#' Check Matrix Name Uniformity in List
#'
#' Verifies if all matrices in a list have identical row and column names.
#' It returns TRUE if uniformity is found, otherwise FALSE.
#'
#' @param matrices List of matrices.
#' @return Logical indicating uniformity of row and column names.
#' @keywords matrix

checkRowColNames <- function(matrices) {
    # Simplify the comparison process
    allEqual <- function(namesList) {
        all(sapply(namesList[-1], function(x) identical(x, namesList[[1]])))
    }

    rowNamesEqual <- allEqual(lapply(matrices, rownames))
    colNamesEqual <- allEqual(lapply(matrices, colnames))

    if ((! rowNamesEqual) || (! colNamesEqual)){
        stop("Rows or columns different among input matrices")
    }
}

#' Parse Metadata From File
#'
#' Reads metadata from a specified file and processes it to handle duplicate
#' rows by aggregating them into a single row based on a unique identifier.
#' The function dynamically identifies the appropriate ID column if not specified.
#' It is designed to be flexible for processing either column (sample) or row (feature) metadata.
#'
#' @param metadata_path Character string specifying the path to the metadata file.
#' @param ids Vector of identifiers (column names or row names) used to match against metadata columns.
#' @param metadata_id_col Optional; character string specifying the column name in the metadata
#'        to be used as the unique identifier. If NULL, the function attempts to
#'        automatically find a suitable column based on `ids`.
#'
#' @return A data frame of processed metadata with duplicate rows aggregated, and row names set to the unique identifier.

parse_metadata <- function(metadata_path, ids, metadata_id_col = NULL){

    metadata <- read_delim_flexible(metadata_path, stringsAsFactors = FALSE, header = TRUE)
    if (is.null(metadata_id_col)){
        metadata_id_col <- findColumnWithAllEntries(ids, metadata)
    }

    # Remove any all-NA columns
    metadata <-  metadata[, colSums(is.na(metadata)) != nrow(metadata)]

    # Allow for duplicate rows by the id column
    metadata <- aggregate(
        . ~ metadata[[metadata_id_col]],
        data = metadata,
        FUN = function(x) paste(unique(x), collapse = ",")
    )[,-1]

    rownames(metadata) <- metadata[[metadata_id_col]]

    metadata[ids,, drop=FALSE]
}

################################################
################################################
## Main script starts here                    ##
################################################
################################################

# Matrices

args_opt <- parse_args('$task.ext.args')
matrix_files <- as.list(strsplit('$matrix_files', ' ')[[1]])

if ('assay_names' %in% names(args_opt)){
    names(matrix_files) <- unlist(strsplit(args_opt[['assay_names']], ',')[[1]])
}else{
    names(matrix_files) <- unlist(lapply(matrix_files, tools::file_path_sans_ext))
}

# Build and verify the main assays list for the summarisedexperiment

assay_list <- lapply(matrix_files, function(m){
    mat <- read_delim_flexible(m, row.names = 1, stringsAsFactors = FALSE)
    mat[,sapply(mat, is.numeric), drop = FALSE]
})

checkRowColNames(assay_list)

# Construct SummarizedExperiment
se <- SummarizedExperiment(
    assays = assay_list
)

# Add column (sample) metadata if provided

if ('$coldata' != ''){
    coldata <- parse_metadata(
        metadata_path = '$coldata',
        ids = colnames(assay_list[[1]]),
        metadata_id_col = args_opt\$coldata_id_col
    )

    colData(se) <- DataFrame(coldata)
}

# Add row (feature) metadata if provided

if ('$rowdata' != ''){
    rowdata <- parse_metadata(
        metadata_path = '$rowdata',
        ids = rownames(assay_list[[1]]),
        metadata_id_col = args_opt\$rowdata_id_col
    )

    rowData(se) <- DataFrame(rowdata)
}

# Write outputs as RDS files

prefix <- tools::file_path_sans_ext(matrix_files[1])
if ('$task.ext.prefix' != 'null'){
    prefix = '$task.ext.prefix'
} else if ('$meta.id' != 'null'){
    prefix = '$meta.id'
}

# Save the SummarizedExperiment object
output_file <- paste0(prefix, ".SummarizedExperiment.rds")
saveRDS(se, file = output_file)

################################################
################################################
## R SESSION INFO                             ##
################################################
################################################

sink(paste(prefix, "R_sessionInfo.log", sep = '.'))
citation("SummarizedExperiment")
print(sessionInfo())
sink()

################################################
################################################
## VERSIONS FILE                              ##
################################################
################################################

r.version <- strsplit(version[['version.string']], ' ')[[1]][3]
summarizedexperiment.version <- as.character(packageVersion('SummarizedExperiment'))

writeLines(
    c(
        '"${task.process}":',
        paste('    bioconductor-summarizedexperiment:', summarizedexperiment.version)
    ),
'versions.yml')

################################################
################################################
################################################
################################################
