% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/cnv_caller_functions_packaged.R
\name{initializeStandards}
\alias{initializeStandards}
\title{initializeStandards}
\usage{
initializeStandards(chromSizeFile, cpgDataFile, cytobandFile)
}
\arguments{
\item{chromSizeFile}{- a tsv chrom sizes file}

\item{cpgDataFile}{- a binned cpg density file (see documentation for how to prepare}

\item{cytobandFile}{- a binned cytoband file for your genome of interest (see documentation for how to prepare)}
}
\description{
This function prepares and preloads the genome data required for CNV analysis.
}
\examples{
initializeStandards("hg38.chrom.sizes","hg38.cpg.1e6.bed","hg38.cytoband.1e6.bed"
}
\keyword{initialisation}
