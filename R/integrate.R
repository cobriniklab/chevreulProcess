#' Split SingleCellExperiment by colData variable
#'
#' @param x SingleCellExperiment object
#' @param f colData variable as a string
#'
#' @return a list of singlecellexperiments name by colData value
#' @export
#'
#' @examples
#' 
#' data(small_example_dataset)
#' splitByCol(small_example_dataset, "batch")
splitByCol <- function(x, f = "batch") {
    f <- colData(x)[[f]]

    i <- split(seq_along(f), f)

    v <- vector(mode = "list", length = length(i))

    names(v) <- names(i)

    for (n in names(i)) {
        v[[n]] <- x[, i[[n]]]
    }

    return(v)
}

#' Merge Small SingleCellExperiment Objects
#'
#' @param ... two or more singlecell objects
#' @param k.filter minimum cell number for integration
#'
#' @return a SingleCellExperiment object
merge_small_sces <- function(..., k.filter = 50) {
    sce_list <- list(...)

    # check if any singlecell objects are too small and if so merge
    # with the first singlecell objects
    sce_dims <- map(sce_list, dim) |>
        map_lgl(~ .x[[2]] < k.filter)

    small_sces <- sce_list[sce_dims]

    sce_list <- sce_list[!sce_dims]

    sce_list[[1]] <- reduce(c(small_sces, sce_list[[1]]), 
                               correctExperiments, PARAM = NoCorrectParam())

    return(sce_list)
}

#' Batch Correct Multiple Single Cell Objects
#'
#' @param sce_list List of two or more SingleCellExperiment objects
#' @param organism human or mouse
#' @param ... extra args passed to sce_reduce_dimensions
#'
#' @return an integrated SingleCellExperiment object
integrate <- function(sce_list, organism = "human", ...) {
    # drop 'colData' fields with same name as 'batchCorrect' output
    sce_list <- map(sce_list, ~ {
        colData(.x)[["batch"]] <- NULL
        return(.x)
    })

    geneCorrected <- correctExperiments(sce_list)
    mainExpName(geneCorrected) <- "integrated"

    geneMerged <- correctExperiments(sce_list, PARAM = NoCorrectParam())
    altExp(geneCorrected, "gene") <- geneMerged

    alt_exp_names <- map(sce_list, altExpNames)

    if (all(map_lgl(alt_exp_names, ~ {
        length(.x) > 0 & "transcript" %in% .x
    }))) {
        # drop 'colData' fields with same name as 'batchCorrect' output
        sce_list <- map(sce_list, ~ {
            colData(.x)[["batch"]] <- NULL
            return(.x)
        })

        transcriptBatches <- map(sce_list, swapAltExp, "transcript")
        transcriptMerged <- correctExperiments(transcriptBatches, 
                                               PARAM = NoCorrectParam())
        altExp(geneCorrected, "transcript") <- transcriptMerged
    }

    geneCorrected <- record_experiment_data(geneCorrected, 
                                            experiment_name = "integrated", 
                                            organism = organism)

    return(geneCorrected)
}

#' Reintegrate (filtered) SingleCellExperiment objects
#'
#' This function takes a SCE object and performs the below steps
#' 1) split by batch
#' 2) integrate
#' 3) run integration pipeline and save
#'
#' @param object A SingleCellExperiment objects
#' @param suffix to be appended to file saved in output dir
#' @param reduction to use default is pca
#' @param ... extra args passed to sce_integrate
#' @export
#'
#' @return a SingleCellExperiment object
reintegrate_sce <- function(object, suffix = "", reduction = "PCA", ...) {
    organism <- metadata(object)$experiment$organism
    experiment_name <- metadata(object)$experiment$experiment_name
    objects <- splitByCol(object, "batch")
    object <- sce_integrate(objects, suffix = suffix, ...)
    object <- record_experiment_data(object, experiment_name, organism)
}
