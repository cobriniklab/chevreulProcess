#' Run Differential Expression
#'
#'
#'
#' @param object a SingleCellExperiment object
#' @param cluster1 cluster 1
#' @param cluster2 cluster 2
#' @param resolution resolution
#' @param diffex_scheme scheme for differential expression
#' @param featureType gene or transcript
#' @param tests t, wilcox, or bimod
#' @export
#' @examples 
#' data("tiny_sce")
#' sce_de(tiny_sce, 
#' colnames(tiny_sce)[1:100], 
#' colnames(tiny_sce)[101:200], 
#' diffex_scheme = "custom")
#'
#' @return a dataframe with differential expression information
sce_de <- function(object, cluster1, cluster2, resolution = 0.2,
                          diffex_scheme = "louvain", featureType = "gene",
                          tests = c("t", "wilcox", "bimod")) {

    data_env <- new.env(parent = emptyenv())
    data("grch38", envir = data_env, package = "chevreulProcess")
    data("grch38_tx2gene", envir = data_env, package = "chevreulProcess")
    grch38 <- data_env[["grch38"]]
    grch38_tx2gene <- data_env[["grch38_tx2gene"]]

    tests <- match.arg(tests)

    if (featureType == "transcript") object <- altExp(object, "transcript")

    if (diffex_scheme == "louvain") {
        if (query_experiment(object, "integrated")) {
            active_experiment <- "integrated"
        } else {
            active_experiment <- "gene"
        }
        colLabels(object) <- colData(object)[[paste0(active_experiment,
                                                     "_snn_res.", resolution)]]
        object <- object[, colLabels(object) %in% c(cluster1, cluster2)]
        colLabels(object) <- factor(colLabels(object))
    } else if (diffex_scheme == "custom") {
        object <- object[, c(cluster1, cluster2)]
        keep_cells <- c(cluster1, cluster2)
        new_idents <- c(rep(1, length(cluster1)), rep(2, length(cluster2)))
        names(new_idents) <- keep_cells
        new_idents <- new_idents[colnames(object)]
        colLabels(object) <- new_idents
        cluster1 <- 1
        cluster2 <- 2
    }
    test_list <- vector("list", length(tests))
    for (test in tests) {
        message(test)
        de <- findMarkers(object, test.type = test)
        if (featureType == "transcript") {
            de_cols <- c("enstxp", "ensgene", "symbol", "p_val" = "p.value",
                         "avg_log2FC", "pct.1", "pct.2", "p_val_adj" = "FDR")
            de <- de[[1]] |>
                as.data.frame() |>
                rownames_to_column("enstxp") |>
                left_join(grch38_tx2gene, by = "enstxp") |>
                left_join(grch38, by = "ensgene")
            if ("summary.logFC" %in% colnames(de)) {
                de <- mutate(de, avg_log2FC = log(exp(summary.logFC), 2))
            }
            de <- select(de, any_of(de_cols))
        } else if (featureType == "gene") {
            de_cols <- c("ensgene", "symbol", 
                         "p_val" = "p.value", "avg_log2FC",
                         "pct.1", "pct.2", "p_val_adj" = "FDR")
            de <- de[[1]] |>
                as.data.frame() |>
                rownames_to_column("symbol") |>
                left_join(grch38, by = "symbol")
            if ("summary.logFC" %in% colnames(de)) {
                de <- mutate(de, avg_log2FC = log(exp(summary.logFC), 2))
            }
            de <- select(de, any_of(de_cols))
        }
        test_list[[match(test, tests)]] <- de
    }
    names(test_list) <- tests
    return(test_list)
}

#' Preprocess Single Cell Object
#'
#' Performs standard pre-processing workflow for scRNA-seq data
#'
#' @param object Assay to use
#' @param scale Perform linear transformation 'Scaling'
#' @param normalize Perform normalization
#' @param features Identify highly variable features
#' @param legacy_settings Use legacy settings
#' @param ... extra args passed to scaling functions
#' @return a preprocessed SingleCellExperiment object
#' @export
#' @examples
#' data("small_example_dataset")
#' sce_preprocess(small_example_dataset)
#'
sce_preprocess <- 
    function(object, scale = TRUE, normalize = TRUE, 
                              features = NULL, legacy_settings = FALSE, ...) {
    clusters <- quickCluster(object)
    object <- computeSumFactors(object, clusters = clusters)
    # summary(sizeFactors(object))

    object <- logNormCounts(object)

    return(object)
}


#' Find All Markers
#'
#' Find all markers at a range of resolutions
#'
#' @param object An object.
#' @param group_by A metadata variable to group by.
#' @param experiment Assay to use, Default "gene".
#' @param ... extra args passed to stash_marker_features
#'
#' @return a SingleCellExperiment object containing marker genes
#' @export 
#' @examples
#' data("small_example_dataset")
#' find_all_markers(small_example_dataset, "gene_snn_res.1")
#'
find_all_markers <- function(object, 
                             group_by = NULL, experiment = "gene", ...) {
    if (is.null(group_by)) {
        meta_cols <- colnames(get_colData(object))
        # resolutions <- meta_cols[grepl(paste0(experiment, "_snn_res."), 
        #                                meta_cols)]
        cluster_index <- grepl(paste0(experiment, "_snn_res."), meta_cols)
        if (!any(cluster_index)) {
            warning("no clusters found in metadata. runnings sce_cluster")
            object <- sce_cluster(object, 
                                     resolution = seq(0.2, 1, by = 0.2))
        }
        clusters <- get_colData(object)[, cluster_index]
        cluster_levels <- map_int(clusters, ~ length(unique(.x)))
        cluster_levels <- cluster_levels[cluster_levels > 1]
        clusters <- select(clusters, one_of(names(cluster_levels)))
        group_by <- names(clusters)
    }
    new_markers <- map(group_by, ~ stash_marker_features(
        object, .x, experiment = experiment, ...))
    names(new_markers) <- group_by
    old_markers <- metadata(object)$markers[!names(
        metadata(object)[["markers"]]) %in% names(new_markers)]
    metadata(object)[["markers"]] <- c(old_markers, new_markers)
    return(object)
}

#' Stash Marker Genes in a SingleCellExperiment Object
#'
#' Marker Genes will be stored in object metadata as `markers`
#'
#' @param group_by A metadata variable to group by
#' @param object A object
#' @param experiment An experiment to use
#' @param top_n Use top n genes, Default 200
#' @param p_val_cutoff p value cut-off, Default value is "0.5"
#'
#' @return a SingleCellExperiment object with marker genes
#' @export
#' @examples
#' data("small_example_dataset")
#' small_example_dataset <- 
#' find_all_markers(small_example_dataset, "gene_snn_res.1")
#' stash_marker_features(small_example_dataset, "gene_snn_res.1")
stash_marker_features <- function(object, group_by, experiment = "gene", 
                                  top_n = 200, p_val_cutoff = 0.5) {
    message("stashing markers for ", group_by)
    markers <- list()
    markers <-
        findMarkers(object, test.type = "t", 
                    groups = colData(object)[[group_by]]) |>
        map(as.data.frame) |>
        map(rownames_to_column, "feature") |>
        bind_rows(.id = "group") |>
        group_by(group) |>
        filter(FDR < p_val_cutoff) |>
        top_n(n = top_n, wt = summary.logFC) |>
        arrange(group, desc(summary.logFC)) |>
        select(Gene.Name = feature, Average.Log.Fold.Change = summary.logFC, 
               Adjusted.pvalue = FDR, Cluster = group) |>
        identity()
    return(markers)
}

#' Run Louvain Clustering at Multiple Resolutions
#'
#' @param object A SingleCellExperiment objects
#' @param resolution Clustering resolution
#' @param custom_clust custom cluster
#' @param reduction Set dimensional reduction object
#' @param algorithm 1
#' @param ... extra args passed to single cell packages
#' @return a SingleCellExperiment object with louvain clusters
#' @export
#' @examples
#' data(small_example_dataset)
#' sce_cluster(small_example_dataset)
#' 
sce_cluster <- function(object = object, resolution = 0.6, 
                           custom_clust = NULL, reduction = "PCA", 
                           algorithm = 1, ...) {
    message(glue("[{format(Sys.time(), '%H:%M:%S')}] Clustering Cells..."))
    if (length(resolution) > 1) {
        for (i in resolution) {
            message(glue("clustering at {i} resolution"))
            cluster_labels <- 
                clusterCells(object,
                             use.dimred = reduction,
                             BLUSPARAM = NNGraphParam(
                                 cluster.fun = "louvain", 
                                 cluster.args = list(resolution = i))
            )
            colData(object)[[glue("gene_snn_res.{i}")]] <- cluster_labels
        }
    } else if (length(resolution) == 1) {
        message(glue("clustering at {resolution} resolution"))
        cluster_labels <- clusterCells(object,
                                       use.dimred = reduction,
                                       BLUSPARAM = NNGraphParam(
                                           cluster.fun = "louvain", 
                                           cluster.args = 
                                               list(resolution = resolution))
        )
        
        
        colData(object)[[glue("gene_snn_res.{resolution}")]] <- cluster_labels
    }
    
    return(object)
}