% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/pipeline.R
\name{sce_integrate}
\alias{sce_integrate}
\title{Run SingleCellExperiment Integration}
\usage{
sce_integrate(
  sce_list,
  resolution = seq(0.2, 1, by = 0.2),
  suffix = "",
  organism = "human",
  annotate_cell_cycle = FALSE,
  annotate_percent_mito = FALSE,
  reduction = "corrected",
  ...
)
}
\arguments{
\item{sce_list}{List of objects to be integrated}

\item{resolution}{Range of resolution}

\item{suffix}{a suffix to be appended to a file save in output dir}

\item{organism}{Default "human"}

\item{annotate_cell_cycle}{whether to score cell cycle phases}

\item{annotate_percent_mito}{logical scalar
whether to annotate mitochondrial percentage}

\item{reduction}{pca, umap, or tsne}

\item{...}{extra args passed to integrate}
}
\value{
an integrated SingleCellExperiment object
}
\description{
Run batch correction, followed by:
\enumerate{
\item stashing of batches in metadata 'batch'
\item clustering with resolution 0.2 to 2.0 in increments of 0.2
\item saving to <proj_dir>/output/sce/\if{html}{\out{<feature>}}\emph{sce}\if{html}{\out{<suffix>}}.rds
}
}
\examples{
data("small_example_dataset")
small_example_dataset |> 
splitByCol("Mutation_Status") |> 
sce_integrate()


}
