
#' Subset by new colData
#'
#' Subset the object using new colData
#'
#' @param colData_path Path to new colData
#' @param object A object
#'
#' @return a SingleCellExperiment object
#'
subset_by_colData <- function(colData_path, object) {
    upload_colData <- read.csv(colData_path, col.names = "sample_id") |>
        filter(!is.na(sample_id) & !sample_id == "sample_id") |>
        mutate(name = sample_id) |>
        column_to_rownames("sample_id") |>
        identity()

    upload_cells <- rownames(upload_colData)

    object <- object[, colnames(object) %in% upload_cells]

    colData(object) <- merge(colData(object), upload_colData, by = "row.names")

    return(object)
}

#' Read in Gene and Transcript SingleCellExperiment Objects
#'
#' @param proj_dir path to project directory
#' @param prefix default "unfiltered"
#'
#' @return a SingleCellExperiment object
load_sce_path <- function(proj_dir = getwd(), prefix = "unfiltered") {
    sce_regex <- paste0(paste0(".*/", prefix, "_sce.rds"))

    sce_path <- path(proj_dir, "output", "singlecellexperiment") |>
        dir_ls(regexp = sce_regex)

    if (!length(sce_path) == 0) {
        return(sce_path)
    }

    stop(sce_path, " does not exist in 
         current working directory ", getwd(), ".",
        call. = FALSE
    )
}


#' Load SingleCellExperiment Files from a single project path
#'
#' @param proj_dir project directory
#' @param ... extra args passed to load_sce_path
#'
#' @return a SingleCellExperiment object
load_sce_from_proj <- function(proj_dir, ...) {
    sce_file <- load_sce_path(proj_dir, ...)
    sce_file <- readRDS(sce_file)
}
