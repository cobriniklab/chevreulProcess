#' Create a database of bigwigfiles
#'
#' Create a sqlite database of bigwig files matching cell ids in objects
#'
#' @param bam_files vector of paths to bam files
#' @param bigwig_db bigwig database
#'
#' @return a path to a bigwig file sqlite database
build_bigwig_db <- function(bam_files, 
                            bigwig_db = "~/.cache/chevreul/bw-files.db") {
    bam_files <- normalizePath(bam_files)

    bigwigfiles <- map_chr(bam_files,
                           ~ bam_to_bigwig(.x,
                                           prefix = path_ext_remove(.x),
                                           overwrite = TRUE))
    
    names(bigwigfiles) <- path_file(bigwigfiles)
    
    bigwigfiles <- 
        bigwigfiles |> 
        enframe("name", "bigWig") |>
        mutate(sample_id = 
                   str_remove(name, "_Aligned.sortedByCoord.out.bw")) |>
        identity()

    con <- dbConnect(SQLite(), dbname = bigwig_db)

    dbWriteTable(con, "bigwigfiles", bigwigfiles, append = TRUE)

    dbDisconnect(con)
}

#' Load Bigwigs
#'
#' Load a tibble of bigwig file paths by cell id
#'
#' @param object A object
#' @param bigwig_db Sqlite database of bigwig files
#'
#' @return a vector of bigwigs file paths
#' @export
load_bigwigs <- function(object, bigwig_db = "~/.cache/chevreul/bw-files.db") {
    con <- dbConnect(SQLite(), dbname = bigwig_db)

    bigwigfiles <- dbReadTable(con, "bigwigfiles") |>
        filter(sample_id %in% colnames(object)) |>
        identity()

    missing_bigwigs <- colnames(object)[!(colnames(object) %in% 
                                              bigwigfiles$sample_id)] |>
        paste(collapse = ", ")

    if(length(missing_bigwigs)>0){
        message(paste0("Sample coverage files ", missing_bigwigs, 
                       "(.bw) do not match samples in object (check file names)"))   
    }

    dbDisconnect(con)

    bigwigfiles <-
        bigwigfiles |>
        filter(sample_id %in% colnames(object))

    return(bigwigfiles)
}

#' Make Bigwig Database
#'
#'
#' @param new_project Project directory
#' @param cache_location Path to cache "~/.cache/chevreul"
#' @param sqlite_db sqlite db containing bw files
#'
#' @return a sqlite database of bigwig files for cells 
#' in a SingleCellExperiment object
make_bigwig_db <- function(new_project = NULL, 
                           cache_location = "~/.cache/chevreul/", 
                           sqlite_db = "bw-files.db") {
    new_bigwigfiles <- dir_ls(new_project, glob = "*.bw", recurse = TRUE)
    
    names(new_bigwigfiles) <- path_file(new_bigwigfiles)
    
    new_bigwigfiles <- 
        new_bigwigfiles |> 
        enframe("name", "bigWig") |>
        mutate(sample_id = 
                   str_remove(name, "_Aligned.sortedByCoord.out.*bw$")) |>
        filter(!str_detect(name, "integrated")) |>
        distinct(sample_id, .keep_all = TRUE) |>
        identity()

    con <- dbConnect(SQLite(), dbname = path(cache_location, sqlite_db))

    all_bigwigfiles <-
        dbReadTable(con, "bigwigfiles") |>
        bind_rows(new_bigwigfiles)

    dbWriteTable(con, "bigwigfiles", all_bigwigfiles, overwrite = TRUE)

    return(all_bigwigfiles)
}
