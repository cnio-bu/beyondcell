#' @title Creates a geneset object
#' @description This function creates a \code{\link[beyondcell]{geneset}}
#' object.
#' @name GenerateGenesets
#' @importFrom qusage read.gmt
#' @importFrom gdata trim
#' @param x A pre-loaded matrix or a path to a GMT file with custom gene sets. 
#' See Details for more information.
#' @param n.genes Number of up and/or down-regulated genes used to compute each
#' signature.
#' @param mode Whether the output \code{geneset} must contain up and/or
#' down-regulated genes. See Details for more information.
#' @param filters If \code{x} is a pre-loaded matrix, you can provide a list of
#' filters to subset it. You can specify which \code{drugs}, sig \code{IDs},
#' mechanisms of action (\code{MoAs}), \code{targets} and/or \code{sources} you
#' are interested in (cap insensitive). You can call
#' \code{\link[beyondcell]{ListFilters}} to check all the available values for
#' these filters. The signatures that pass \strong{ANY} of them are included in
#' the output.
#' @param include.pathways Logical. Return \code{beyoncell}'s pre-computed
#' signatures for functional pathways?
#' @details \code{x} can be:
#' \itemize{
#' \item{A pre-loaded matrix:} {Either \code{PSc}, \code{SSc} or \code{DSS}.}
#' \item{A path to a GMT file:} {A file that contains custom gene sets. Each
#' gene set must have an "_UP" or "_DOWN" suffix.}
#' }
#' In addition, \code{mode} can be:
#' \itemize{
#' \item{\code{"up"}:} {To compute the signatures using only up-regulated
#' genes.}
#' \item{\code{"down"}:} {To compute the signatures using only down-regulated
#' genes.}
#' \item{\code{c("up", "down")} :} {To compute the signatures using both up and
#' down-regulated genes.}
#' }
#' If \code{x} is a path to a GMT file, \code{mode} is deprecated and the names
#' of all gene sets must end in "_UP" or "_DOWN" to indicate the \code{mode} of
#' each one.
#' @return A \code{geneset} object.
#' @examples
#' @export

GenerateGenesets <- function(x, n.genes = 250, mode = c("up", "down"),
                             include.pathways = TRUE) {
  # --- Global Checks ---
  # Check if x is a path to a GMT file.
  message('Reading gmt file...')
  gmt.file <- tryCatch(qusage::read.gmt(x), error = function(cond) {
    message(paste0('The GMT file does not seem to exist: ', x, ' .'))
    stop(paste('x must be either a pre-loaded matrix (PSc, SSc or DSS), a
                ranked matrix or the path to a GMT file.'))
  }, warning = function(cond) {
    message(paste('Warning when reading the GMT file. Here is the',
                  'original warning message:'))
    warning(paste0(cond))
    return(suppressWarnings(qusage::read.gmt(x)))
  })
  ### Check for duplicated gene sets.
  upper.gmt.names <- toupper(names(gmt.file))
  if (anyDuplicated(upper.gmt.names) != 0) {
    duplicated.gene.sets <- unique(names(gmt.file)[duplicated(upper.gmt.names)])
    stop(paste0('The GMT file contains duplicated gene set\'s: ',
                paste0(duplicated.gene.sets, collapse = ", "), '.'))
  }
  
  # Check n.genes and mode.
  if (any(!(mode %in% c("up", "down")))) stop('Incorrect mode.')
  mode <- sort(unique(mode), decreasing = TRUE)

  ### Number of genes.
  if (!identical(n.genes, 250)) warning('x is a GMT file, n.genes is deprecated.')
  ### Mode in GMT files.
  n.up <- length(unique(grep(pattern = "_UP$", x = upper.gmt.names)))
  n.down <- length(unique(grep(pattern = "_DOWN$", x = upper.gmt.names)))
  if (n.up + n.down != length(names(gmt.file))) {
    stop('All gene sets\' names in the GMT file must end in "_UP" or "_DOWN".')
  } else {
    if (n.up > 0 & n.down > 0) {
      if (!identical(mode, c("up", "down"))) {
        mode <- c("up", "down")
        warning(paste('The GMT file includes UP and DOWN gene sets. mode',
                        'changed to c("up", "down").'))
      }
    } else if (n.up > 0) {
      if (mode != "up") {
        mode <- "up"
        warning('The GMT file only includes UP gene sets. mode changed to "up".')
      }
    } else if (n.down > 0) {
      if (mode != "down") {
        mode <- "down"
        warning('The GMT file only includes DOWN gene sets. mode changed to "down".')
      }
    }
  }
  # Check include.pathways.
  if (length(include.pathways) != 1 | !is.logical(include.pathways)) {
    stop('include.pathways must be TRUE or FALSE.')
  }
  
  # --- Code ---
  ### Genes.
  unique.gene.sets <- unique(gsub(pattern = "_UP$|_DOWN$", replacement = "",
                                  x = names(gmt.file), ignore.case = TRUE))
  genes <- setNames(lapply(unique.gene.sets, function(sig) {
    l <- list()
    if (toupper(paste0(sig, "_UP")) %in% upper.gmt.names) {
      l <- c(l, list(up = gmt.file[[match(toupper(paste0(sig, "_UP")),
                                          table = upper.gmt.names)]]))
    }
    if (toupper(paste0(sig, "_DOWN")) %in% upper.gmt.names) {
      l <- c(l, list(down = gmt.file[[match(toupper(paste0(sig, "_DOWN")),
                                            table = upper.gmt.names)]]))
    }
    return(l)
  }), unique.gene.sets)
  # Drug IDs.
  info <- data.frame()
  # Pathways.
  if (include.pathways) {
    paths <- lapply(pathways, function(p) p[names(p)[mode %in% names(p)]])
  } else {
    paths <- list()
  }
  # Output.
  return(geneset(genelist = c(genes, paths), n.genes = n.genes,
                 mode = mode, info = info, inverse.score = inverse.score))
}

#' @title Returns all the possible values for the specified filter
#' @description This function returns all the available values for \code{drugs},
#' \code{IDs}, \code{MoAs}, \code{targets} or \code{sources} filters in
#' \code{\link[beyondcell]{GenerateGenesets}} function.
#' @name ListFilters
#' @param entry Either \code{"drugs"}, \code{"IDs"}, \code{"MoAs"},
#' \code{"targets"} or \code{"sources"}.
#' @return All the possible values for the specified \code{entry}.
#' @export

ListFilters <- function(entry) {
  # --- Checks and Code ---
  if (entry == "drugs") {
    out <- sort(unique(drugInfo$drugs), decreasing = FALSE)
  } else if (entry == "IDs") {
    out <- sort(unique(drugInfo$IDs), decreasing = FALSE)
  } else if (entry == "MoAs") {
    out <- sort(unique(drugInfo$MoAs), decreasing = FALSE)
  } else if (entry == "targets") {
    out <- sort(unique(drugInfo$targets), decreasing = FALSE)
  } else if (entry == "sources") {
    out <- sort(unique(drugInfo$sources), decreasing = FALSE)
  } else {
    stop("Incorrect entry.")
  }
  return(out)
}
