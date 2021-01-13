#' An S4 class to represent the genesets of each drug
#'
#' @slot genelist A list of drug and functional signatures (\code{pathways})
#' with up and/or down-regulated genes.
#' @slot n.genes Argument passed to
#' \code{\link[GenerateGenesets]{GenerateGenesets}}.
#' @slot mode Argument passed to
#' \code{\link[GenerateGenesets]{GenerateGenesets}}.
#' @slot info If \code{GenerateGenesets}' input is a pre-loaded matrix,
#' \code{info} will be a \code{data.frame} with the correspondences between
#' \code{sig_ids} and drug names, as well as the source of the pre-loaded
#' matrices (LINCS, CTRP, GDSC or CCLE). If \code{GenerateGenesets}' input is a
#' path to a GMT file, \code{info} will be empty.
#' @slot comparison Argument passed to
#' \code{\link[GenerateGenesets]{GenerateGenesets}}.

geneset <- setClass("geneset", slots = list(genelist = "list",
                                            n.genes = "numeric",
                                            mode = "character",
                                            info = "data.frame",
                                            comparison = "character"))

#' An S4 class to represent the beyondcell scores for each cell and signature.
#'
#' @slot scaled (Subsetted and/or regressed) scaled beyondcell scores.
#' @slot normalized (Subsetted and/or regressed) normalized beyondcell scores.
#' @slot data Original normalized beyondcell scores, without subsetting or
#' regression.
#' @slot switch.point (Subsetted and/or regressed) scaled beyondcell score for
#' which the normalized score in \code{@@data} is 0 (one switch point per
#' signature).
#' @slot ranks List of dataframes with the statistics (switch point,
#' mean, median, sd, variance, min, max, proportion of \code{NaN}, residuals'
#' mean and ranks of each signature in the beyondcell object. This slot is
#' filled using the function \code{\link[bcRanks]{bcRanks}}.
#' @slot expr.matrix Expression matrix used to compute the beyondcell scores.
#' @slot meta.data Contains information about each cell (including the
#' therapeutic clusters and \code{Seurat}'s \code{@@metadata}).
#' @SeuratInfo List with information about the input \code{Seurat} object,
#' including the \code{@@reductions}.
#' @slot background (Subsetted and/or regressed) normalized beyondcell scores
#' obtained using DSS signatures. Useful to compute the UMAP reduction and the
#' therapeutic clusters when the number of drug signatures is low.
#' @slot reductions A list of dimensional reductions for this object.
#' @slot regression A list with the order of subset and regression steps
#' performed on the beyondcell object and the variables used for regression.
#' @slot n.genes Argument passed to
#' \code{\link[GenerateGenesets]{GenerateGenesets}}.
#' @slot mode Argument passed to
#' \code{\link[GenerateGenesets]{GenerateGenesets}}.
#' @slot thres Argument \code{expr.thres} passed to
#' \code{\link[bcScore]{bcScore}}.

beyondcell <- setClass("beyondcell",
                       slots = list(scaled = "matrix", normalized = "matrix",
                                    data = "matrix", switch.point = "numeric",
                                    ranks = "list", expr.matrix = "matrix",
                                    meta.data = "data.frame", SeuratInfo = "list",
                                    background = "matrix", reductions = "list",
                                    regression = "list", n.genes = "numeric",
                                    mode = "character", thres = "numeric"))
