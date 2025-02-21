% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/cnv_caller_functions_packaged.R
\name{filterCells}
\alias{filterCells}
\title{filterCells}
\usage{
filterCells(inputMatrix, minimumSegments = 40, minDensity = 0, signalSDcut = 2)
}
\arguments{
\item{inputMatrix}{- collapsed, normalised matrix}

\item{minimumSegments}{the minimum segments required to keep a cell}

\item{minDensity}{minimum segment density to consider useful (defaults to zero)}

\item{signalSDcut}{Number of standard deviations average signal of a cell must be from an average per-cell signal to be considered an outlier (to identify putative doublets; default: 2)}
}
\description{
This function filters the collapsed chromosome signal density matrix based on predefined parameters (minimum non-zero chromosome segments)
}
\examples{
filterCells(inputMatrix,minimumSegments=41)
}
\keyword{CNV}
\keyword{filter}
