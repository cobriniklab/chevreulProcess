#' Get cell metadata from a given object
#' 
#' Get cell metadata
#'
#' @param object a SingleCellExperiment object
#'
#' @return dataframe containing object metadata
#' @export
#' @examples
#' 
#' data(small_example_dataset)
#' get_colData(small_example_dataset)
get_colData <- function(object) {
    colData(object) |>
        as.data.frame()
}

#' Set cell metadata
#' 
#' Set cell metadata from a given object
#'
#' @param object a SingleCellExperiment object
#' @param meta a dataframe containing object metadata
#'
#' @return a SingleCellExperiment object with new colData
#' @export
#' @examples
#' 
#' data(small_example_dataset)
#' new_meta <- data.frame(row.names = colnames(small_example_dataset))
#' new_meta$example <- "example"
#' set_colData(small_example_dataset, new_meta)
set_colData <- function(object, meta) {
    colData(object) <- DataFrame(meta)
    return(object)
}

#' Get object metadata
#'
#' @param object a SingleCellExperiment object
#'
#' @return variable features from a SingleCellExperiment object
#' @importFrom S4Vectors metadata
get_sce_metadata <- function(object) {
    metadata(object)
}

#' Get variable features
#'
#' @param object a SingleCellExperiment object
#' @param experiment "gene" or "transcript"
#'
#' @return variable features from a SingleCellExperiment object
#' @export
#' @examples
#' 
#' data(small_example_dataset)
#' get_variable_features(small_example_dataset)
get_variable_features <- function(object, experiment = "gene") {
    if (experiment == mainExpName(object)) {
        getTopHVGs(object)
    } else if (experiment %in% altExpNames(object)) {
        getTopHVGs(altExp(object, experiment))
    }
}

#' Get feature names
#'
#' @param object a SingleCellExperiment object
#' @param experiment "gene" or "transcript"
#'
#' @return variable features from a SingleCellExperiment object
#' @export
#' @examples
#' 
#' data(small_example_dataset)
#' get_features(small_example_dataset)
get_features <- function(object, experiment = "gene") {
    if (experiment == mainExpName(object)) {
        rownames(object)
    } else if (experiment %in% altExpNames(object)) {
        rownames(altExp(object, experiment))
    }
}

#' Get Feature Types
#'
#' @param object a SingleCellExperiment object
#'
#' @return vector of feature types in an object
#' @export
#'
#' @examples
#' 
#' data(small_example_dataset)
#' get_feature_types(small_example_dataset)
get_feature_types <- function(object) {
    sort(c(mainExpName(object), altExpNames(object)))
}

#' Set Feature Types
#'
#' @param object a SingleCellExperiment object
#' @param feature_type feature type
#' @return a SingleCellExperiment object with assigned feature type
#' @export
#'
#' @examples
#' 
#' data(small_example_dataset)
#' set_feature_type(small_example_dataset, "transcript")
set_feature_type <- function(object, feature_type) {
    if (feature_type %in% altExpNames(object)) {
        object <- swapAltExp(object, feature_type, saved = mainExpName(object), withColData = TRUE)
    }
    return(object)
}

#' Retrieve Assay
#'
#' @param object a SingleCellExperiment object
#' @param experiment an experiment name
#'
#' @return Main or alt experiment in a SingleCellExperiment object
retrieve_experiment <- function(object, experiment) {
    if (experiment %in% mainExpName(object)) {
        return(object)
    } else if (experiment %in% altExpNames(object)) {
        return(altExp(object, experiment))
    }
}

#' Query Experiment
#'
#' @param object a SingleCellExperiment object
#' @param experiment an experiment name
#'
#' @return logical scalar indicating if experiment is present in object
#' @export
#' @examples 
#' data(small_example_dataset)
#' query_experiment(small_example_dataset, "gene")
query_experiment <- function(object, experiment) {
    return(experiment %in% c(mainExpName(object), altExpNames(object)))
}
