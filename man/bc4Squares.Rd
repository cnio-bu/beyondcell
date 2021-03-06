% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/Visualization.R
\name{bc4Squares}
\alias{bc4Squares}
\title{Drawns a 4 squares plot}
\usage{
bc4Squares(
  bc,
  idents,
  lvl = NULL,
  top = 3,
  topnames = NULL,
  force = 1,
  alpha = 0.7,
  pt.size = 3,
  ...
)
}
\arguments{
\item{bc}{\code{beyondcell} object.}

\item{idents}{Name of the metadata column of interest.}

\item{lvl}{Character vector with the \code{idents}' level(s) of interest. If
\code{lvl = NULL}, all levels will be plotted.}

\item{top}{Number of top drugs per quadrant to be labelled.}

\item{topnames}{Character vector with additional interesting drugs or 
pathways to be labeled (either their names or sig IDs).}

\item{force}{(From \code{\link[ggrepel]{ggrepel}}) Force of repulsion between
overlapping text labels. Defaults to 1.}

\item{alpha}{Transparency level between 1 (not transparent) and 0 (fully
transparent).}

\item{pt.size}{Point size.}

\item{...}{Other arguments passed to \code{\link[ggrepel]{geom_text_repel}}.}
}
\value{
A list of \code{ggplot2} objects, one per \code{lvl}.
}
\description{
This function drawns a 4 square plot of the drug signatures
present in a \code{\link[beyondcell]{beyondcell}} object. A 4 squares plot
consists in a scatter plot of the residuals' means (x axis) vs the switch
points (y axis) of a specified cluster (either a therapeutic cluster or a
group defined by experimental condition or phenotype). 4 quadrants are
highlighted: the top-left and bottom-right corners contain the drugs to which
all selected cells are least/most sensistive, respectively. The centre
quadrants show the drugs to which these cells are differentially insensitive
or sensitive when compared to the other clusters.

x cut-offs: first and last deciles; y cut-offs: 0.1, 0.4, 0.6 and 0.9.
}
\details{
This function returns a list of \code{\link[ggplot2]{ggplot2}}
objects, one per each \code{lvl}. Note that residuals' means are different
for each level while swicth points are signature-specific. So, x axis will
vary and y axis will remain constant accross all plots.
}
