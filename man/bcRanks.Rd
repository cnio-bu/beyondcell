% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/Ranks.R
\name{bcRanks}
\alias{bcRanks}
\title{Ranks the signatures from most sensitive to least sensitive}
\usage{
bcRanks(bc, idents = NULL, extended = TRUE)
}
\arguments{
\item{bc}{\code{\link[beyondcell]{beyondcell}} object.}

\item{idents}{Name of the metadata column of interest. If
\code{idents = NULL}, the function computes the ranks using all cells. If
\code{idents != NULL}, the signatures' ranks are computed for each level in
\code{idents}.}

\item{extended}{If \code{extended = TRUE}, this function returns the switch
point, mean, median, standard deviation, variance, min, max, proportion of
\code{NaN} and residuals' mean per signature. If \code{extended = FALSE},
this function returns only the switch point, mean and residuals' mean.}
}
\value{
A \code{beyondcell} object with the results in a new entry of
\code{bc@ranks}: \code{bc@ranks[["general"]]} (if \code{idents = NULL}) or
\code{bc@ranks[[idents]]} (if \code{idents != NULL}).
}
\description{
This function computes the beyondcell score's (BCS) statistics
of each signature and ranks them according to the switch point and mean.
}
