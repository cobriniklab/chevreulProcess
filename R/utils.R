#' Get Transcripts in object
#'
#' Get transcript ids in objects for one or more gene of interest
#'
#' @param object A SingleCellExperiment object
#' @param gene Gene of interest
#'
#' @return transcripts constituting a 
#' gene of interest in a SingleCellExperiment object
get_transcripts_from_sce <- function(object, gene) {
    transcripts <- genes_to_transcripts(gene)

    transcripts <- transcripts[transcripts %in% get_features(object, 
                                                             "transcript")]
}

#' Record Experiment Metadata
#'
#' Records miscellaneous data
#' @param object A object
#' @param experiment_name name of the experiment
#' @param organism human or mouse
#'
#' @return a SingleCellExperiment object
#' @export
#' @examples
#' 
#' data(small_example_dataset)
#' record_experiment_data(small_example_dataset)
#'
record_experiment_data <- function(object, 
                                   experiment_name = "default_experiment", 
                                   organism = "human") {
    if (!requireNamespace("SingleCellExperiment", quietly = TRUE)) {
        stop("Package 'SingleCellExperiment' needed for this function to work. 
             Please install it.",
            call. = FALSE
        )
    }

    organism <- metadata(object)[["experiment"]][["organism"]] %||% organism

    experiment_name <- 
        metadata(object)[["experiment"]][["experiment_name"]] %||% experiment_name

    message(glue("[{format(Sys.time(), '%H:%M:%S')} 
                 Logging Technical Details..."))
    experiment <- list(
        experiment_name = experiment_name,
        organism = organism
    )
    experiment$date_of_export <- Sys.Date()
    experiment$date_of_analysis <- Sys.Date()

    experiment$parameters <- list(
        gene_nomenclature = "gene_symbol",
        discard_genes_expressed_in_fewer_cells_than = 10,
        keep_mitochondrial_genes = TRUE,
        variables_to_regress_out = "nCount_RNA",
        number_PCs = 30,
        tSNE_perplexity = 30,
        cluster_resolution = seq(0.2, 1, by = 0.2)
    )
    experiment$filtering <- list(
        UMI_min = 50,
        genes_min = 10
    )
    experiment$sessionInfo <- list(
        capture.output(sessionInfo())
    )

    if (!is.null(objectVersion(object))) {
        experiment$SingleCellExperiment_version <- objectVersion(object)
    }

    experiment$chevreul_version <- packageVersion("chevreulProcess")

    metadata(object)[["experiment"]] <- experiment

    return(object)
}

#' Calculate Read Count Metrics for a object
#'
#' Recalculate counts/features per cell for a object
#'
#' @param object A SingleCellExperiment object
#'
#' @return a SingleCellExperiment object with nfeatures and 
#' ngenes stored in metadata
#' @export
#' @examples
#' 
#' data(small_example_dataset)
#' sce_calcn(small_example_dataset)
sce_calcn <- function(object) {
    object <- addPerCellQC(object)
    main_exp_name <- mainExpName(object) %||% "main"
    object[[glue("nFeature_{main_exp_name}")]] <- object$detected
    object[[glue("nCount_{main_exp_name}")]] <- object$sum

    for (alt_exp_name in altExpNames(object)) {
        altExp(object, alt_exp_name) <- addPerCellQC(altExp(object, 
                                                            alt_exp_name))
        altExp(object, alt_exp_name)[[glue("nFeature_{alt_exp_name}")]] <- 
            altExp(object, alt_exp_name)[["detected"]]
        altExp(object, alt_exp_name)[[glue("nCount_{alt_exp_name}")]] <- 
            altExp(object, alt_exp_name)[["sum"]]
    }

    return(object)
}

#' Propagate Metadata Changes
#'
#' @param meta updated metadata
#' @param object a SingleCellExperiment object
#'
#' @return a SingleCellExperiment object
#' @export
#' @examples
#'
#' 
#' data(small_example_dataset)
#' new_meta <- data.frame(row.names = colnames(small_example_dataset))
#' new_meta$example <- "example"
#'
#' propagate_spreadsheet_changes(new_meta, small_example_dataset)
propagate_spreadsheet_changes <- function(meta, object) {
    meta <- meta |>
        rownames_to_column("cell") |>
        mutate(across(contains("snn"), as.factor)) |>
        mutate(across(where(is.ordered), ~ as.factor(as.character(.x)))) |>
        column_to_rownames("cell") |>
        identity()

    colData(object) <- DataFrame(meta)

    return(object)
}

#' Create a database of chevreul projects
#'
#' Create a database containing chevreul projects
#'
#' @param cache_location Path to cache "~/.cache/chevreul"
#' @param sqlite_db Database to be created
#' @param verbose print messages
#'
#' @return a sqlite database with SingleCellExperiment objects
create_project_db <- function(cache_location = "~/.cache/chevreul", 
                              sqlite_db = "single-cell-projects.db", 
                              verbose = TRUE) {
    if (!dir.exists(cache_location)) {
        dir.create(cache_location)
    }
    con <- dbConnect(SQLite(), path(cache_location, sqlite_db))
    projects_tbl <- tibble(project_name = character(), 
                           project_path = character(), 
                           project_slug = character(), 
                           project_type = character(), )
    message(glue(
        "building table of chevreul projects at {path(cache_location, 
        sqlite_db)}"))
    tryCatch({
        dbWriteTable(con, "projects_tbl", projects_tbl)
    }, warning = function(w) {
        message(sprintf("Warning in %s: %s", deparse(w[["call"]]), 
                        w[["message"]]))
    }, error = function(e) {
        message("projects db already exists!")
    }, finally = {
    })
    dbDisconnect(con)
}

#' Update a database of chevreul projects
#'
#' Add new/update existing projects to the database by recursing fully
#'
#' @param projects_dir The project directory to be updated
#' @param cache_location Path to cache "~/.cache/chevreul"
#' @param sqlite_db sqlite db
#' @param verbose print messages
#'
#' @return a sqlite database with SingleCellExperiment objects
update_project_db <- function(
        projects_dir = NULL,
        cache_location = "~/.cache/chevreul",
        sqlite_db = "single-cell-projects.db",
        verbose = TRUE) {
    if (!dir.exists(cache_location)) {
        dir.create(cache_location)
    }

    con <- dbConnect(SQLite(), path(cache_location, sqlite_db))

    projects_tbl <-
        dir_ls(projects_dir, glob = "*.here", recurse = TRUE, 
               fail = FALSE, all = TRUE) |>
        path_dir()
    
    names(projects_tbl) <- path_file(projects_tbl)
    
    projects_tbl <- 
        projects_tbl |>
        enframe("project_name", "project_path") |>
        mutate(project_slug = str_remove(project_name, "_proj$")) |>
        mutate(project_type = path_file(path_dir(project_path))) |>
        identity()

    current_projects_tbl <-
        dbReadTable(con, "projects_tbl") |>
        filter(file.exists(project_path)) |>
        filter(!project_path %in% projects_tbl$project_path) |>
        bind_rows(projects_tbl) |>
        distinct(project_path, .keep_all = TRUE)

    dbWriteTable(con, "projects_tbl", projects_tbl, overwrite = TRUE)

    dbDisconnect(con)
}

#' Update a database of chevreul projects
#'
#' Append projects to database
#'
#' @param new_project_path new project path
#' @param cache_location Path to cache "~/.cache/chevreul"
#' @param sqlite_db sqlite db
#' @param verbose print messages
#' @return a sqlite database with SingleCellExperiment objects
append_to_project_db <- function(
        new_project_path,
        cache_location = "~/.cache/chevreul",
        sqlite_db = "single-cell-projects.db",
        verbose = TRUE) {
    if (!dir.exists(cache_location)) {
        dir.create(cache_location)
    }

    con <- dbConnect(SQLite(), path(cache_location, sqlite_db))

    projects_tbl <-
        new_project_path
    
    names(projects_tbl) <- path_file(projects_tbl)
    
    projects_tbl <- 
        projects_tbl |> 
        enframe("project_name", "project_path") |>
        mutate(project_slug = str_remove(project_name, "_proj$")) |>
        mutate(project_type = path_file(path_dir(project_path))) |>
        identity()

    current_projects_tbl <-
        dbReadTable(con, "projects_tbl") |>
        filter(file.exists(project_path)) |>
        filter(!project_path %in% projects_tbl$project_path) |>
        bind_rows(projects_tbl) |>
        distinct(project_path, .keep_all = TRUE)

    dbWriteTable(con, "projects_tbl", current_projects_tbl, overwrite = TRUE)

    dbDisconnect(con)
}

#' Read a database of chevreul projects
#'
#' Reads database of chevreul projects to a data frame
#'
#' @param cache_location Path to cache "~/.cache/chevreul"
#' @param sqlite_db sqlite db
#' @param verbose print messages
#'
#' @return a tibble with SingleCellExperiment objects
#
read_project_db <- function(
        cache_location = "~/.cache/chevreul",
        sqlite_db = "single-cell-projects.db",
        verbose = TRUE) {
    if (!dir.exists(cache_location)) {
        dir.create(cache_location)
    }

    con <- dbConnect(SQLite(), path(cache_location, sqlite_db))

    current_projects_tbl <-
        dbReadTable(con, "projects_tbl")

    dbDisconnect(con)

    return(current_projects_tbl)
}

#' Retrieve Metadata from Batch
#'
#' @param batch batch
#' @param projects_dir path to project dir
#' @param db_path path to .db file
#'
#' @return a tibble with cell level metadata from a SingleCellExperiment object
metadata_from_batch <- function(
        batch, projects_dir = "/dataVolume/storage/single_cell_projects",
        db_path = "single-cell-projects.db") {
    mydb <- dbConnect(SQLite(), path(projects_dir, db_path))

    projects_tbl <- dbReadTable(mydb, "projects_tbl") |>
        filter(!project_type %in% c("integrated_projects", "resources"))

    dbDisconnect(mydb)

    metadata <-
        projects_tbl |>
        filter(project_slug == batch) |>
        pull(project_path) |>
        path("data") |>
        dir_ls(glob = "*.csv") |>
        identity()
}


#' Save object to <project>/output/sce/<feature>_sce.rds
#'
#' @param object a SingleCellExperiment object
#' @param prefix a prefix for saving
#' @param proj_dir path to a project directory
#'
#' @return a path to an rds file containing a SingleCellExperiment object
#'
#'
#'
save_sce <- function(object, prefix = "unfiltered", proj_dir = getwd()) {
    sce_dir <- path(proj_dir, "output", "singlecellexperiment")

    dir.create(sce_dir)

    sce_path <- path(sce_dir, paste0(prefix, "_sce.rds"))

    message(glue("saving to {sce_path}"))
    saveRDS(object, sce_path)

    return(object)
}
