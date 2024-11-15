#' Regress SingleCellExperiment Object by Given Set of Genes
#'
#' @param object A object
#' @export
#'
#' @return a SingleCellExperiment object with features regressed
regress_cell_cycle <- function(object) {
    message("regressing objects by cell cycle")

    if (dim(object)[2] < 100) {
        ctrl <- dim(object)[2] / 10
    }
    if (!"Phase" %in% colnames(colData(object))) {
        object <- annotate_cell_cycle(object)
    }
    if (!any(str_detect(names(assays(object)), pattern = ".*corrected*"))) {
        dec.nocycle <- modelGeneVar(object, block = colData(object)[["Phase"]])
        reg.nocycle <- regressBatches(object, batch = colData(object)[["Phase"]])
        
        assay(object, "corrected") <- assay(reg.nocycle, "corrected")
        
        object <- runPCA(object,
            exprs_values = "corrected",
            subset_row = getTopHVGs(dec.nocycle, prop = 0.1)
        )

        original_experiment <- mainExpName(object)

        res_names <- grepl(glue("{original_experiment}_snn_res."), 
                           colnames(colData(object)))
        resolutions <- 
            str_extract(colnames(colData(object))[res_names], "[0-9].*$")
        object <- runTSNE(x = object, dimred = "PCA", n_dimred = seq(30))
        object <- runUMAP(x = object, dimred = "PCA", n_dimred = seq(30))
        object <- sce_cluster(object, resolution = resolutions)
    }
    return(object)
}





#' Annotate Cell Cycle
#'
#' Annotate Cell Cycle for Gene and Transcript SingleCellExperiment Objects
#'
#' @param object A SingleCellExperiment object
#'
#' @return a SingleCellExperiment object
annotate_cell_cycle <- function(object) {
    
    data_env <- new.env(parent = emptyenv())
    data("cc.genes.cyclone", envir = data_env, package = "chevreulProcess")
    cc.genes.cyclone <- data_env[["cc.genes.cyclone"]]
    
    assignments <- cyclone(object, cc.genes.cyclone, 
                           gene.names = rownames(object))
    colData(object)[colnames(assignments$scores)] <- assignments$scores
    colData(object)["Phase"] <- assignments$phases
    return(object)
}