#' Run SingleCellExperiment Integration
#'
#' Run batch correction, followed by:
#' 1) stashing of batches in metadata 'batch'
#' 2) clustering with resolution 0.2 to 2.0 in increments of 0.2
#' 3) saving to <proj_dir>/output/sce/<feature>_sce_<suffix>.rds
#'
#' @param suffix a suffix to be appended to a file save in output dir
#' @param sce_list List of objects to be integrated
#' @param resolution Range of resolution
#' @param organism Default "human"
#' @param annotate_cell_cycle whether to score cell cycle phases
#' @param reduction pca, umap, or tsne
#' @param ... extra args passed to integrate
#' @param annotate_percent_mito logical scalar
#' whether to annotate mitochondrial percentage
#'
#' @return an integrated SingleCellExperiment object
sce_integrate <- function(sce_list, resolution = seq(0.2, 2.0, by = 0.2), suffix = "", organism = "human", annotate_cell_cycle = FALSE, annotate_percent_mito = FALSE, reduction = "PCA", ...) {
    experiment_names <- names(sce_list)

    organisms <- case_when(
        grepl("Hs", experiment_names) ~ "human",
        grepl("Mm", experiment_names) ~ "mouse"
    )

    names(organisms) <- experiment_names

    organisms[is.na(organisms)] <- organism

    integrated_sce <- integrate(sce_list, organism = organism, ...)

    integrated_sce <- sce_reduce_dimensions(integrated_sce, ...)

    # cluster merged objects
    integrated_sce <- sce_cluster(integrated_sce, resolution = resolution, algorithm = algorithm, reduction = reduction, ...)

    experiment <- "gene"
    integrated_sce <- find_all_markers(integrated_sce, experiment = experiment)

    #   enriched_sce <- tryCatch(getEnrichedPathways(integrated_sce), error = function(e) e)
    #   enrichr_available <- !any(class(enriched_sce) == "error")
    #   if(enrichr_available){
    #     integrated_sce <- enriched_sce
    #   }

    # annotate cell cycle scoring to objects
    if (annotate_cell_cycle) {
        integrated_sce <- annotate_cell_cycle(integrated_sce, ...)
    }

    # annotate mitochondrial percentage in object metadata
    if (annotate_percent_mito) {
        integrated_sce <- add_percent_mito(integrated_sce, ...)
    }

    # annotate excluded cells
    # integrated_sce <- annotate_excluded(integrated_sce, excluded_cells)

    return(integrated_sce)
}

# sce_process ------------------------------

#' Run SingleCellExperiment Pipeline
#'
#' This functions allows you to preprocess, cluster and reduce dimensions 
#' for one SingleCellExperiment object.
#'
#' @param object A SingleCellExperiment object
#' @param experiment Assay of interest in SingleCellExperiment object
#' @param resolution Resolution for clustering cells. Default set to 0.6.
#' @param reduction Dimensional reduction object
#' @param organism Organism
#' @param ... extra parameters passed to internal functions
#'
#' @return a processed SingleCellExperiment object
#' @export
#' @examples
#' data(small_example_dataset)
#' sce_process(small_example_dataset)
#' 
sce_process <- function(object, experiment = "gene", resolution = 0.6, reduction = "PCA", organism = "human", ...) {
    object <- sce_preprocess(object, scale = TRUE, ...)
    for (experiment in altExpNames(object)) {
        altExp(object, experiment) <- sce_preprocess(altExp(object, experiment), scale = TRUE, ...)
    }

    # PCA
    object <- sce_reduce_dimensions(object, ...)

    object <- sce_cluster(object = object, resolution = resolution, reduction = reduction, ...)

    object <- find_all_markers(object, experiment = "gene")

    # if (feature == "gene"){
    #   enriched_sce <- tryCatch(getEnrichedPathways(object), error = function(e) e)
    #   enrichr_available <- !any(class(enriched_sce) == "error")
    #   if(enrichr_available){
    #     object <- enriched_sce
    #   }
    # }

    # annotate cell cycle scoring to objects
    object <- annotate_cell_cycle(object, ...)

    # annotate mitochondrial percentage in object metadata
    object <- add_percent_mito(object, ...)

    return(object)
}

#' Dimensional Reduction
#'
#' Run PCA, TSNE and UMAP on a singlecell objects
#' perplexity should not be bigger than 3 * perplexity < nrow(X) - 1, 
#' see details for interpretation
#'
#' @param object A SingleCellExperiment object
#' @param experiment Experiment of interest to be processed
#' @param ... Extra parameters passed to sce_reduce_dimensions
#'
#' @return a SingleCellExperiment object with embeddings
sce_reduce_dimensions <- function(object, experiment = "gene", ...) {

    num_samples <- dim(object)[[2]]
    if (num_samples < 50) {
        npcs <- num_samples - 1
    } else {
        npcs <- 50
    }
    if ("gene" == experiment) {
        object <- runPCA(x = object, subset_row = getTopHVGs(stats = object), 
                         ncomponents = npcs, ...)
    } else {
        object <- runPCA(x = object, altexp = experiment, 
                         subset_row = getTopHVGs(stats = object), 
                         ncomponents = npcs, ...)
    }
    
    if ((ncol(object) - 1) > 3 * 30) {
        if ("gene" == experiment) {
            object <- runTSNE(x = object, dimred = "PCA", n_dimred = seq(30))
        } else {
            object <- runTSNE(x = object, altexp = experiment, dimred = "PCA", 
                              n_dimred = seq(30))
        }
        if ("gene" == experiment) {
            object <- runUMAP(x = object, dimred = "PCA", n_dimred = seq(30))
        } else {
            object <- runUMAP(x = object, altexp = experiment, dimred = "PCA", 
                              n_dimred = seq(30))
        }
    }

    return(object)
}

#' Give a new project name to a SingleCellExperiment object
#'
#' @param object A SingleCellExperiment object
#' @param new_name New name to assign
#'
#' @return a renamed SingleCellExperiment object
#' @export
#' @examples
#' 
#' 
#' data(small_example_dataset)
#' rename_sce(small_example_dataset, "new_name")
rename_sce <- function(object, new_name) {
    metadata(object)["project.name"] <- new_name
    return(object)
}
