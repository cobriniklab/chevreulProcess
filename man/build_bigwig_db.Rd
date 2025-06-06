% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/wiggleplotr.R
\name{build_bigwig_db}
\alias{build_bigwig_db}
\title{Create a database of bigwigfiles}
\usage{
build_bigwig_db(bam_files, bigwig_db = "~/.cache/chevreul/bw-files.db")
}
\arguments{
\item{bam_files}{vector of paths to bam files}

\item{bigwig_db}{bigwig database}
}
\value{
a path to a bigwig file sqlite database
}
\description{
Create a sqlite database of bigwig files matching cell ids in objects
}
