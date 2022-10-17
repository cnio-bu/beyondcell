#' @title Creates a geneset object
#' @description This function creates a \code{\link[beyondcell]{geneset}}
#' object.
#' @name GenerateGenesets
#' @importFrom qusage read.gmt
#' @importFrom gdata trim
#' @param x A path to a GMT file with custom gene sets. 
#' See Details for more information.
#' @param perform.reversal Logical. Return \code{beyoncell}'s pre-computed
#' signatures as if the signature was reversed.
#' @details \code{x} can be:
#' \itemize{
#' \item{A path to a GMT file:} {A file that contains custom gene sets. Each
#' gene set must have an "_UP" or "_DOWN/_DN" suffix.}
#' }
#' @return A \code{geneset} object.
#' @examples
#' @export

GenerateGenesets <- function(x, perform.reversal = FALSE) {
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

  ### Mode in GMT files.
  n.up <- length(unique(grep(pattern = "_UP$", x = upper.gmt.names)))
  n.down <- length(unique(grep(pattern = "_DOWN$|_DN$", x = upper.gmt.names)))
  if (n.up + n.down != length(names(gmt.file))) {
    stop('All gene sets\' names in the GMT file must end in "_UP" or "_DOWN/_DN".')
  } 

    # Check perform.reversal.
  if (length(perform.reversal) != 1 | !is.logical(perform.reversal)) {
    stop('perform.reversal must be TRUE or FALSE.')
  }
  # --- Code ---
  ### Genes.
  unique.gene.sets <- unique(gsub(pattern = "_UP$|_DOWN$|_DN$", replacement = "",
                                  x = names(gmt.file), ignore.case = TRUE))
  genes <- setNames(lapply(unique.gene.sets, function(sig) {
    l <- list()
    if (toupper(paste0(sig, "_UP")) %in% upper.gmt.names) {
      l <- c(l, list(up = gmt.file[[match(toupper(paste0(sig, "_UP")),
                                          table = upper.gmt.names)]]))
    }
    if (toupper(paste0(sig, "_DOWN")) %in% upper.gmt.names |
        toupper(paste0(sig, "_DN")) %in% upper.gmt.names ){
      # Match will return NA if no matches available.
      sig_dn_down <- paste0(sig, c("_DOWN", "_DN"))
      sig_dn_down <- toupper(sig_dn_down)
      match_id <- match(sig_dn_down, table = upper.gmt.names)
      valid_id <- match_id[!is.na(match_id)]

      l <- c(l, list(down = gmt.file[[valid_id]]))
    }
    return(l)
  }), unique.gene.sets)

  # Output.
  return(geneset(genelist = genes,
  n.genes = NaN, # User defined genesets have no n.genes
  mode = "None",  # User defined genesets have no mode
  info = data.frame(), # User defined genesets have an empty info. slot.
  inverse.score = perform.reversal))
}

#' @title Load beyondcell collection
#' @description This function creates a \code{\link[beyondcell]{geneset}}
#' object from one of the default collections available.
#' @name GetCollection
#' @importFrom gdata trim
#' @param x Either PSc, SSc or DSS.
#' See Details for more information.
#' @param n.genes Number of up and/or down-regulated genes used to compute each
#' signature.
#' @param mode Whether the output \code{geneset} must contain up and/or
#' down-regulated genes. See Details for more information.
#' @param filters If \code{x} is a pre-loaded collection, you can provide a list of
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
#' @return A \code{geneset} object.
#' @examples
#' @export

GetCollection <- function(x, n.genes = 250, mode = c("up", "down"),
                          filters = list(drugs = NULL, IDs = NULL, MoAs = NULL,
                                         targets = NULL, studies = NULL),
                          include.pathways = TRUE){

   # --- Global Checks ---
  # Check if x is a preloaded collection
  collection_name <- deparse(substitute(x))
  is.D <- c(identical(x, PSc), identical(x, SSc), identical(x, DSS))
  n.max <- 500 # Our colls. have at max. 500 genes

  if (!any(is.D)) {
   stop(paste('x must be either PSc, SSc or DSS.'))
  } 

  # Check n.genes and mode.
  if (any(!(mode %in% c("up", "down")))) stop('Incorrect mode.')
  mode <- sort(unique(mode), decreasing = TRUE)

  ### Number of genes.
  if (length(n.genes) != 1 | n.genes[1]%%1 != 0 | n.genes[1] < 1) {
    stop('n.genes must be a positive integer.')
    } else if (n.genes > n.max) {
      stop(paste0('n.genes exceeds the maximum number of genes in signature (',
                  n.max, ').'))
    }

  # Check filters.
  filters.names <- c("drugs", "IDs", "MoAs", "targets", "studies")
  selected.filters <- names(filters)
  if (any(!(selected.filters %in% filters.names))) stop('Invalid names in filters.')
  filters_class <- sapply(filters, is.null) | sapply(filters, is.character)
    if (any(!filters_class)) {
      stop(paste0('Incorrect value for filter\'s entry: "',
                  paste0(selected.filters[!filters_class], collapse = ", "),
                  '". You must provide a character vector.'))
    }
    selected.filters <- selected.filters[!sapply(filters, is.null)]

    # Check include.pathways.
  if (length(include.pathways) != 1 | !is.logical(include.pathways)) {
    stop('include.pathways must be TRUE or FALSE.')
  }

   # --- Code ---
  # Subset pre-loaded collections...
  if (collection_name == "DSS"){
    dss_sigs <- names(DSS@genelist)
    info <- subset(drugInfo[["IDs"]],
                   subset = drugInfo[["IDs"]]$IDs %in% dss_sigs)
                   } else {
    info <- subset(drugInfo[["IDs"]], 
                   subset = drugInfo[["IDs"]]$collections == collection_name
                 )}
                 
  inverse.score <- FALSE
  if (identical(x, PSc) | identical(x, DSS)) {
    inverse.score <- TRUE # When using PSc/DDS, inverse the sign of the BCS.
  }
  ### Filters.
  if (length(selected.filters) == 0) {
    ids <- unique(info$IDs)
    } else {
      ids <- unique(unlist(lapply(selected.filters, function(y) {
        tryCatch(suppressWarnings(GetIDs(values = filters[[y]], filter = y)),
                 error = function(cond) character())
      })))
      warnings <- unlist(lapply(selected.filters, function(z) {
        tryCatch(GetIDs(values = filters[[z]], filter = z),
                 error = function(cond) {
                   err <- paste0(z, ": ", paste0(filters[[z]],
                                                 collapse = ", "), ".\n")
                   return(err)
                 }, warning = function(cond) {
                   warn <- as.character(cond)
                   warn.values <- strsplit(sapply(strsplit(warn, split = ": "),
                                                  `[[`, 3), split = ", ")
                   return(paste0(z, ": ", warn.values))
                 })
      }))
      warnings <- warnings[!startsWith(warnings, prefix = "sig_")]
      if (length(ids) == 0) {
        stop('Couldn\'t find signatures that matched any of the filters.')
      } else if (length(warnings) > 0) {
        warning(paste('The following filters\' values yielded no results:\n',
                      paste0("   - ", warnings, " ", collapse = "")))
      }
    }
  ### Genes.
  genes <- lapply(ids, function(sig) {
    l <- list(up = x@genelist[[sig]]$up[1:n.genes],
              down = x@genelist[[sig]]$down[1: n.genes]
              )
    return(l)
  })
  names(genes) <- ids
   # Drug IDs.
  info <- subset(info, subset = info$IDs %in% ids)
  info <- aggregate(.~ IDs, data = info, na.action = NULL, FUN = function(rw) {
    paste(na.omit(unique(rw)), collapse = ", ")
  })
  info <- info[order(info$IDs, decreasing = FALSE), ]
  
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
    out <- sort(unique(drugInfo[["Synonyms"]]$drugs), decreasing = FALSE)
  } else if (entry == "IDs") {
    out <- sort(unique(drugInfo[["IDs"]]$IDs), decreasing = FALSE)
  } else if (entry == "MoAs") {
    out <- sort(unique(drugInfo[["MoAs"]]$MoAs), decreasing = FALSE)
  } else if (entry == "targets") {
    out <- sort(unique(drugInfo[["Targets"]]$targets), decreasing = FALSE)
  } else if (entry == "sources") {
    out <- sort(unique(drugInfo[["IDs"]]$studies), decreasing = FALSE)
  } else {
    stop("Incorrect entry.")
  }
  return(out)
}

