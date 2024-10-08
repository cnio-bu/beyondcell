% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/Basics.R
\name{FindDrugs}
\alias{FindDrugs}
\title{Returns a dataframe with information about the input drugs}
\usage{
FindDrugs(bc, x, na.rm = TRUE)
}
\arguments{
\item{bc}{\code{\link[beyondcell]{beyondcell}} object.}

\item{x}{A character vector with drug names and/or sig IDs.}

\item{na.rm}{Logical. Should \code{x} entries with no available drug 
information be removed from the final output?}
}
\value{
A \code{data.frame}.
}
\description{
This function searches the input drugs in the pre-loaded
\code{beyondcell} matrices and returns a dataframe with drug information,
including drug synonyms and MoAs.
}
\details{
The output \code{data.frame} has the following columns:
\itemize{
\item{\code{original.names}}: Input drug names.
\item{\code{bc.names}}: Drug names used in \code{bc}.
\item{\code{preferred.drug.names}}: Standard drug names.
\item{\code{drugs}}: Other drug names.
\item{\code{IDs}}: Signature IDs.
\item{\code{preferred.and.sigs}}: \code{preferred.drug.names} (or alternatively
\code{bc.names}) and \code{IDs}. Used as title in \code{beyondcell} plots.
\item{\code{MoAs}}: Mechanism(s) of action.
}
}
