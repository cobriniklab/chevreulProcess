% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/pipeline.R
\name{sce_process}
\alias{sce_process}
\title{Run SingleCellExperiment Pipeline}
\usage{
sce_process(
  object,
  experiment = "gene",
  resolution = 0.6,
  reduction = "PCA",
  organism = "human",
  process = TRUE,
  ...
)
}
\arguments{
\item{object}{A SingleCellExperiment object}

\item{experiment}{Assay of interest in SingleCellExperiment object}

\item{resolution}{Resolution for clustering cells. Default set to 0.6.}

\item{reduction}{Dimensional reduction object}

\item{organism}{Organism}

\item{process}{whether to run dimensional reduction and clustering}

\item{...}{extra parameters passed to internal functions}
}
\value{
a processed SingleCellExperiment object
}
\description{
This functions allows you to preprocess, cluster and reduce dimensions
for one SingleCellExperiment object.
}
\examples{
data(tiny_sce)
sce_process(tiny_sce, process = FALSE)

}
