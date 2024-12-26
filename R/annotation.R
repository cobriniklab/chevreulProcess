#' Gene Symbols to Ensembl Transcript Ids
#'
#' convert hgnc gene symbols to ensembl transcript ids
#'
#' @param symbols character vector of gene symbols
#'
#' @return a vector of transcripts
#' @export
#' @importFrom GenomicFeatures transcripts
#'
#' @examples
#'
#' genes_to_transcripts("NRL")
genes_to_transcripts <- function(symbols) {
    
        feature_table <- transcripts(EnsDb.Hsapiens.v86,
            columns = c("gene_name", "gene_biotype", "gene_id"),
            return.type = "DataFrame")

    feature_table[(feature_table[["gene_name"]] %in% symbols), "tx_id"]
}

#' Ensembl Transcript Ids to Gene Symbols
#'
#' Convert ensembl transcript ids to hgnc gene symbols
#'
#' @param transcripts human transcripts
#'
#' @return a vector of gene symbols
#' @export
#'
#' @examples
#'
#' NRL_transcripts_hs <-
#' c("ENST00000359842", "ENST00000470566", "ENST00000465764")
#'
#' transcripts_to_genes(transcripts = NRL_transcripts_hs)
#'
transcripts_to_genes <- function(transcripts) {
    
    data_env <- new.env(parent = emptyenv())
    data("grch38", envir = data_env, package = "chevreulProcess")
    data("grch38_tx2gene", envir = data_env, package = "chevreulProcess")
    grch38 <- data_env[["grch38"]]
    grch38_tx2gene <- data_env[["grch38_tx2gene"]]

    tibble(enstxp = transcripts) |>
        left_join(grch38_tx2gene, by = "enstxp") |>
        left_join(grch38, by = "ensgene") |>
        pull("symbol") |>
        identity()
}

#' Annotate percent mitochondrial reads per cell
#'
#'  Add a Percentage of Mitochondrial Read Count Categorical Variable to the
#'  Object (based on nCount_RNA)
#'
#' @param object A object
#' @param experiment gene
#'
#' @return a single cell object with
#' cell metadata column containing mitochondrial percentage
#'
add_percent_mito <- function(object, experiment = "gene") {
    is.mito <- grepl("^MT-*", rownames(object))
    object <- addPerCellQCMetrics(object, subsets = list(Mito = is.mito))
    return(object)
}
