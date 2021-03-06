% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/Format.R
\name{center_scale_colour_stepsn}
\alias{center_scale_colour_stepsn}
\title{Creates a centred sequential binned colour gradient}
\usage{
center_scale_colour_stepsn(
  x,
  colorscale,
  alpha = 0.7,
  na.value = "grey50",
  limits = c(NA, NA),
  center = NULL,
  breaks = 0.1
)
}
\arguments{
\item{x}{A numeric vector. It can contain \code{NA}s.}

\item{colorscale}{A vector with 5 colours that can be obtained using
\code{\link[beyondcell]{get_colour_steps}}.}

\item{alpha}{Transparency level between 1 (not transparent) and 0 (fully
transparent).}

\item{na.value}{Colour to use for missing values.}

\item{limits}{Vector with the desired limits.}

\item{center}{A single number indicating the centre of the \code{colorscale}.
If \code{center = NULL} (default), the centre is set to the middle point of
\code{x}.}

\item{breaks}{A single number indicating the break size of the
\code{colorscale}. Alternatively, it can be a vector with the desired breaks
(which don't have to be symmetric or equally distributed).}
}
\value{
A centred sequential binned colour gradient that can be used to
colour \code{\link[ggplot2]{ggplot2}} objects.
}
\description{
This function creates a sequential binned colour gradient
(low-mid-high) centred around \code{center}.
}
