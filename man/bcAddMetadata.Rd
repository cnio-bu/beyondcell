% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/Manipulation.R
\name{bcAddMetadata}
\alias{bcAddMetadata}
\title{Add new metadata to an existing beyondcell object}
\usage{
bcAddMetadata(bc, metadata)
}
\arguments{
\item{bc}{\code{beyondcell} object.}

\item{metadata}{Matrix or dataframe with metadata to add. Rownames should be
cell names and colnames should not be already present in
\code{bc@meta.data}.}
}
\value{
Returns a \code{beyondcell} object with updated metadata.
}
\description{
This function adds new metadata to an existing
\code{\link[beyondcell]{beyondcell}} object.
}
