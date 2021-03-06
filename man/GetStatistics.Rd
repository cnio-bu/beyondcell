% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/Basics.R
\name{GetStatistics}
\alias{GetStatistics}
\title{Computes the BCS' statistics and ranks}
\usage{
GetStatistics(bc, signatures, cells, pb, total, i, n.rows, extended)
}
\arguments{
\item{bc}{\code{\link[beyondcell]{beyondcell}} object.}

\item{signatures}{Vector with the names of the signatures of interest.}

\item{cells}{Vector with the names of the cells of interest.}

\item{pb}{\code{\link[utils]{txtProgressBar}}.}

\item{total}{Number of iterations to complete the \code{pb}.}

\item{i}{Iteration number. Used to increase the \code{pb}.}

\item{n.rows}{Number of signatures. Used to increase the \code{pb}.}

\item{extended}{If \code{extended = TRUE}, this function returns the switch
point, mean, median, sd, variance, min, max, proportion of \code{NaN} and
residuals' mean per signature. If \code{extended = FALSE}, this function
returns only the switch point, mean and residuals' mean.}
}
\value{
A \code{data.frame} with the BCS' statistics and ranks.
}
\description{
This function computes the beyondcell scores' (BCS) statistics
and ranks returned by \code{\link[beyondcell]{bcRanks}}.
}
