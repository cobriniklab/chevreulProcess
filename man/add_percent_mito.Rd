% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/annotation.R
\name{add_percent_mito}
\alias{add_percent_mito}
\title{Annotate percent mitochondrial reads per cell}
\usage{
add_percent_mito(object, experiment = "gene")
}
\arguments{
\item{object}{A object}

\item{experiment}{gene}
}
\value{
a single cell object with
cell metadata column containing mitochondrial percentage
}
\description{
Add a Percentage of Mitochondrial Read Count Categorical Variable to the
Object (based on nCount_RNA)
}
