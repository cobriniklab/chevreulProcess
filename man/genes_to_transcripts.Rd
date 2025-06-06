% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/annotation.R
\name{genes_to_transcripts}
\alias{genes_to_transcripts}
\title{Gene Symbols to Ensembl Transcript Ids}
\usage{
genes_to_transcripts(symbols)
}
\arguments{
\item{symbols}{character vector of gene symbols}
}
\value{
a vector of transcripts
}
\description{
convert hgnc gene symbols to ensembl transcript ids
}
\examples{

genes_to_transcripts("NRL")
}
