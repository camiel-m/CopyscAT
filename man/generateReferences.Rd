% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/cnv_caller_functions_packaged.R
\name{generateReferences}
\alias{generateReferences}
\title{generateReferences
Generate reference files from UCSC for genome of interest
Returns CPG, chromosome size and cytoband data into target directory}
\usage{
generateReferences(
  genomeObject,
  genomeText = "hg38",
  tileWidth = 1e+06,
  outputDir = "~"
)
}
\arguments{
\item{genomeObject}{BSgenome object of interest}

\item{genomeText}{Shorthand for genome in UCSC (defaults to hg38)}

\item{tileWidth}{Width of tiles - (default: 1e6)}

\item{outputDir}{Output directory to save reference files to (default: ~)}
}
\description{
generateReferences
Generate reference files from UCSC for genome of interest
Returns CPG, chromosome size and cytoband data into target directory
}
\keyword{CNV}
\keyword{reference}
