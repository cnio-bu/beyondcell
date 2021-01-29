#' @title Creates a geneset object
#' @description This function generates up and/or down genelists for each drug
#' and pathway.
#' @name GenerateGenesets
#' @importFrom qusage read.gmt
#' @importFrom gdata trim
#' @param x A pre-loaded matrix, a numeric matrix or a path to a GMT file with
#' custom genesets.
#' @param n.genes Number of top and/or down genes in the desired genelists.
#' @param mode \code{"up"}, \code{"down"} or \code{c("up", "down")}. See Details
#' for more information.
#' @param filters If \code{x} is a pre-loaded matrix, you can provide a list of
#' filters to subset these matrices. You can specify which drug names, sig IDs,
#' mechanisms of action, target genes and sources you are interested in (cap
#' insensitive). You can call \code{\link[ListFilters]{ListFilters}} function to
#' check all the available values for these filters. The signatures that pass
#' \strong{ANY} of these filters are included in the output.
#' @param comparison \code{"treated_vs_control"} or
#' \code{"control_vs_treated"}. See Details for more information.
#' @param include.pathways \code{TRUE} (default) or \code{FALSE}. Whether or
#' not return \code{beyoncell}'s pre-computed genesets for functional pathways.
#' @details \code{x} can be:
#' \itemize{
#' \item{A pre-loaded matrix:} {Either \code{\link[PSc]{PSc}},
#' \code{\link[SSc]{SSc}} or \code{\link[DSS]{DSS}}.}
#' \item{A numeric matrix:} {A matrix with genes as rows and signatures as
#' columns that contains some type of numerical value such a t-stat or a LFC to
#' rank the genes accordingly.}
#' \item{A path to a GMT file:} {A file that contains custom genesets. Each
#' geneset must have an "_UP" or "_DOWN" suffix.}
#' }
#' In addition, \code{mode} can be:
#' \itemize{
#' \item{\code{"up"}:} {To return over-expressed genes in the signatures.}
#' \item{\code{"down"}:} {To return under-expressed genes in the signatures.}
#' }
#' If \code{x} is a path to a GMT file, \code{mode} is deprecated and the names
#' of all genesets must end in "_UP" or "_DOWN" to indicate the mode of each
#' one.
#'
#' Finally, \code{comparison} can be:
#' \itemize{
#' \item{\code{"treated_vs_control"}:} {(\code{PSc} and \code{DSS} like) When
#' the numeric values or the genesets in the GMT file were obtained from a
#' comparison between drug treated and untreated cells.}
#' \item{\code{"control_vs_treated"}:} {(\code{SSc} like) When the numeric
#' values or the genesets in the GMT file were obtained from a comparison
#' between drug sensitive and resistant cells.}
#' }
#' When \code{x} is a pre-loaded matrix, \code{comparison} is set automatically.
#' @return A \code{\link[geneset]{geneset}} object with up and/or down genes for
#' each drug and pathway.
#' @examples
#' @export

GenerateGenesets <- function(x, n.genes = 250, mode = c("up", "down"),
                             filters = list(drugs = NULL, IDs = NULL, MoA = NULL,
                                            targets = NULL, source = NULL),
                             comparison = NULL, include.pathways = TRUE) {
  # --- Global Checks ---
  # Check if x is a pre-loaded matrix, a numeric matrix or a path to a GMT file.
  is.D <- c(identical(x, PSc), identical(x, SSc), identical(x, DSS))
  if (any(is.D)) {
    type <- "pre-loaded matrix"
    message(paste('Reading', c("PSc", "SSc", "DSS")[is.D], 'signatures...'))
    n.max <- 500
  } else if (is.matrix(x) | is.data.frame(x)) {
    type <- "matrix"
    message('Reading input matrix...')
    ### Check if x is numeric.
    if (!is.numeric(x)) stop('x must be a numeric matrix.')
    x <- as.matrix(x)
    ### Check if there are missing values.
    if (sum(is.na(x)) > 0) {
      warning('x contains NAs that will not be used to compute the genesets.')
    }
    ### Check if there are duplicated values in the same column.
    dup.values <- apply(x, 2, FUN = function(y) any(duplicated(na.omit(y))))
    if (is.null(colnames(x))) colnames(x) <- 1:ncol(x)
    if(is.null(rownames(x))) rownames(x) <- 1:nrow(x)
    if (any(dup.values)) {
      warning(paste0('The following columns contain duplicated values:',
                     paste0(colnames(x)[dup.values], collapse = ", "), '.'))

    }
    n.max <- max(apply(x, 2, FUN = function(z) length(na.omit(z))))
  } else {
    type <- "gmt"
    message('Reading gmt file...')
    gmt.file <- tryCatch(qusage::read.gmt(x), error = function(cond) {
      message(paste0('The GMT file does not seem to exist: ', x, ' .'))
      stop(paste('x must be either a pre-loaded matrix (PSc, SSc or DSS) or',
                 'the path to a GMT file.'))
    }, warning = function(cond) {
      message(paste('Warning when reading the GMT file. Here is the',
                    'original warning message:'))
      warning(paste0(cond))
      return(suppressWarnings(qusage::read.gmt(x)))
    })
    ### Check for duplicated genesets.
    if (anyDuplicated(names(gmt.file)) != 0) {
      duplicated_genesets <- unique(names(gmt.file)[duplicated(names(gmt.file))])
      stop(paste0('The GMT file contains duplicated geneset\'s names: ',
                  paste0(duplicated_genesets, collapse = ", "), '.'))
    }
  }
  # Check n.genes and mode.
  if (!all(mode %in% c("up", "down"))) stop('Incorrect mode.')
  mode <- sort(unique(mode), decreasing = TRUE)
  if (type != "gmt") {
    ### Number of genes.
    if (length(n.genes) != 1 | !is.numeric(n.genes)) {
      stop('n.genes must be a single number.')
    }
    if (n.genes%%1 != 0 | n.genes < 1) {
      stop('n.genes must be a positive integer.')
    } else if (n.genes > n.max) {
      stop(paste0('n.genes exceeds the maximum number of genes in signature (',
                  n.max, ').'))
    }
  } else {
    ### Number of genes.
    if (!identical(n.genes, 250)) warning('x is a GMT file, n.genes is deprecated.')
    ### Mode in GMT files.
    n.up <- length(unique(grep("_UP$", names(gmt.file), ignore.case = TRUE)))
    n.down <- length(unique(grep("_DOWN$", names(gmt.file), ignore.case = TRUE)))
    if (n.up + n.down != length(names(gmt.file))) {
      stop('All geneset in the GMT file names must be finished in _UP or _DOWN.')
    } else {
      if (n.up > 0 & n.down > 0) {
        if (!identical(mode, c("up", "down"))) {
          mode <- c("up", "down")
          warning(paste('The GMT file includes UP and DOWN genesets. mode',
                        'changed to c("up", "down").'))
        }
      } else if (n.up > 0) {
        if (mode != "up") {
          mode <- "up"
          warning('The GMT file only includes UP genesets. mode changed to "up".')
        }
      } else if (n.down > 0) {
        if (mode != "down") {
          mode <- "down"
          warning('The GMT file only includes DOWN genesets. mode changed to "down".')
        }
      }
    }
  }
  # Check filters and comparison.
  filters_names <- c("drugs", "IDs", "MoA", "targets", "source")
  selected_filters <- names(filters)
  if (!all(selected_filters %in% filters_names)) stop('Invalid names in filters.')
  if (type != "pre-loaded matrix") {
    ### Filters.
    filters_class <- sapply(filters, is.null)
    if (any(!filters_class)) {
      warning('x is not a pre-loaded matrix, filters is deprecated.')
    }
    ### Comparison.
    if (is.null(comparison)) {
      stop(paste('Comparison must be either "treated_vs_control" or',
                 '"control_vs_treated".'))
    } else if (length(comparison) != 1 |
               !(comparison[1] %in% c("treated_vs_control", "control_vs_treated"))) {
      stop(paste('Comparison must be either "treated_vs_control" or',
                 '"control_vs_treated".'))
    }
  } else {
    ### Filters.
    filters_class <- sapply(filters, is.null) | sapply(filters, is.character)
    if (any(!filters_class)) {
      stop(paste0('Incorrect value for filter\'s entry: "',
                  paste0(selected_filters[!filters_class], collapse = ", "),
                  '". You must provide a character vector.'))
    }
    selected_filters <- selected_filters[!sapply(filters, is.null)]
    ### Comparison.
    if (is.null(comparison)) {
      if(is.D[2]) comparison <- "control_vs_treated"
      else comparison <- "treated_vs_control"
    } else {
      if (length(comparison) != 1 |
          !(comparison[1] %in% c("treated_vs_control", "control_vs_treated"))) {
        stop('Incorrect comparison.')
      }
      if (is.D[2] & comparison != "control_vs_treated") {
        comparison <- "control_vs_treated"
        warning('x = SSc, comparison changed to "control_vs_treated".')
      } else if (!is.D[2] & comparison != "treated_vs_control") {
        comparison <- "treated_vs_control"
        warning(paste0('x = ', c("PSc", "SSc", "DSS")[is.D], ', comparison ',
                       'changed to "treated_vs_control".'))
      }
    }
  }
  # Check include.pathways.
  if (length(include.pathways) != 1 | !is.logical(include.pathways)) {
    stop('include.pathways must be TRUE or FALSE.')
  }
  # --- Code ---
  # If x is a pre-loaded matrix...
  if (type == "pre-loaded matrix") {
    ### sig IDs.
    if (is.D[1]) {
      info <- subset(drugInfo, drugInfo$Source == "LINCS")
    } else if (is.D[2]) {
      info <- subset(drugInfo, drugInfo$Source != "LINCS")
    } else if (is.D[3]) {
      info <- subset(drugInfo, drugInfo$sig_id %in% DSS[[1]]$sig_id)
      x <- PSc # DSS is a subset of PSc
    }
    if (length(selected_filters) == 0) {
      ids <- unique(info$sig_id)
    } else {
      ids <- c()
      warn <- c() # Warnings
      ### Filters.
      if("drugs" %in% selected_filters) {
        out <- FilteredIDS(info, "Name", gdata::trim(filters$drugs), "drugs")
        ids <- c(ids, out[[1]])
        warn <- c(warn, out[[2]])
      }
      if("IDs" %in% selected_filters) {
        out <- FilteredIDS(info, "sig_id", gdata::trim(filters$IDs), "IDs")
        ids <- c(ids, out[[1]])
        warn <- c(warn, out[[2]])
      }
      if ("MoA" %in% selected_filters) {
        out <- FilteredIDS(info, "MoA", gdata::trim(filters$MoA), "MoAs")
        ids <- c(ids, out[[1]])
        warn <- c(warn, out[[2]])
      }
      if ("targets" %in% selected_filters) {
        out <- FilteredIDS(info, "Target", gdata::trim(filters$targets),
                           "target genes")
        ids <- c(ids, out[[1]])
        warn <- c(warn, out[[2]])
      }
      if ("source" %in% selected_filters) {
        out <- FilteredIDS(info, "Source", gdata::trim(filters$source),
                           "source databases")
        ids <- c(ids, out[[1]])
        warn <- c(warn, out[[2]])
      }
      ids <- unique(ids)
      if (length(ids) == 0) {
        stop('Couldn\'t find drugs that matched any of the filters.')
      } else {
        if (!all(is.null(warn))) sapply(warn[!is.null(warn)], function(w) warning(w))
      }
    }
    genes <- lapply(ids, function(sig) {
      l <- list(up = x[["up"]][1:n.genes, sig], down = x[["down"]][1:n.genes, sig])
      return(l)
    })
    names(genes) <- ids
    # Else if x is a numeric matrix...
  } else if (type == "matrix") {
    genes <- apply(x, 2, FUN = function(sig) {
      l <- list()
      if ("up" %in% mode) {
        up <- na.omit(rownames(x)[order(sig, decreasing = TRUE,
                                        na.last = NA)[1:n.genes]])
        l <- c(l, list(up = up))
      }
      if ("down" %in% mode) {
        down <- na.omit(rownames(x)[order(sig, na.last = NA)[1:n.genes]])
        l <- c(l, list(down = down))
      }
      return(l)
    })
    # Else if x is a GMT file.
  } else if (type == "gmt") {
    unique_genesets <- unique(gsub("_UP$|_DOWN$", "", names(gmt.file)))
    genes <- setNames(lapply(unique_genesets, function(sig) {
      entry <- gmt.file[grep(paste0("^", sig), names(gmt.file), value = TRUE)]
      return(Genelist(entry, mode = mode))
    }), unique_genesets)
  }
  # Drug IDs.
  if (type == "pre-loaded matrix") {
    info <- subset(info, info$sig_id %in% ids)
    info <- aggregate(.~ sig_id, data = info, na.action = NULL, FUN = function(rw) {
      paste(na.omit(unique(rw)), collapse = ", ")
    })
    info <- info[order(info$sig_id), ]
  } else {
    info <- data.frame()
  }
  # Pathways.
  if (include.pathways) {
    paths <- lapply(pathways, function(p) p[names(p)[mode %in% names(p)]])
  } else {
    paths <- list()
  }
  # Output.
  return(geneset(genelist = c(genes, paths), n.genes = n.genes,
                 mode = mode, info = info, comparison = comparison))
}

#' @title Returns all the possible values for the specified filter
#' @description This function returns all the available values for \code{drugs},
#' \code{sig_ids}, \code{MoAs}, \code{targets} or \code{source} filters in
#' \code{\link[GenerateGenesets]{GenerateGenesets}} function.
#' @name ListFilters
#' @param entry Either \code{"drugs"}, \code{"IDs"}, \code{"MoA"},
#' \code{"targets"} or \code{"source"}.
#' @return All the possible values for the specified \code{entry}.
#' @export

ListFilters <- function(entry) {
  # --- Checks and Code ---
  if (entry == "drugs") {
    out <- sort(unique(drugInfo$Name))
  } else if (entry == "IDs") {
    out <- sort(unique(drugInfo$sig_id))
  } else if (entry == "MoA") {
    out <- sort(unique(drugInfo$MoA))
  } else if (entry == "targets") {
    out <- sort(unique(drugInfo$Target))
  } else if (entry == "source") {
    out <- sort(unique(drugInfo$Source))
  } else {
    stop("Incorrect entry.")
  }
  return(out)
}

#' @title Returns the \code{sig_ids} that match the specified filters
#' @description This function is meant to be used inside
#' \code{\link[GenerateGenesets]{GenerateGenesets}}. It subsets \code{infodf}
#' to select only the entries that match the specified \code{filter} and returns
#' the corresponding \code{sig_ids}.
#' @name FilteredIDS
#' @param infodf \code{data.frame} with all drug information.
#' @param col Column name to subset by.
#' @param filter User-supplied filtering vector for either drugs, MoA, target
#' genes or source database.
#' @param filtername String to be printed with the warning (in case some of the
#' \code{filter} elements are not found in \code{infodf}).
#' @return A vector with the \code{sig_ids} that match the \code{filter}
#' elements.
#' @export

FilteredIDS <- function(infodf, col, filter, filtername) {
  # --- Checks ---
  # Check infodf.
  if (class(infodf) != "data.frame") {
    stop('infodf must be a data.frame')
  }
  if (!("sig_id" %in% colnames(infodf))) {
    stop('infodf must contain a "sig_id" column.')
  }
  # Check col.
  if (length(col) != 1) {
    stop('You must specify a single col.')
  }
  if (is.character(col) & !(col %in% colnames(infodf))) {
    stop(paste('col =', col, 'is not a column of infodf.'))
  }
  if (is.numeric(col) & (col < 1 | col > ncol(infodf))) {
    stop(paste('col = ', col, 'is out of range.'))
  }
  # Check filter.
  if (length(filter) < 1 | !is.character(filter)) {
    stop('filter must be a character vector.')
  }
  # Check filtername.
  if (length(filtername) != 1 | !is.character(filtername)) {
    stop('filtername must be a string.')
  }
  # --- Code ---
  upper.filter <- toupper(filter)
  selected <- subset(infodf, toupper(infodf[[col]]) %in% upper.filter)
  if (col == "Name") {
    synonyms <- subset(infodf, toupper(infodf[["Preferred_Name"]]) %in%
                         unique(toupper(selected[["Preferred_Name"]])))
    selected <- unique(rbind(selected, synonyms))
  }
  ids <- unique(selected$sig_id)
  not_found <- filter[!(upper.filter %in% toupper(infodf[[col]]))]
  if (length(not_found) > 0) {
    w <- paste(length(not_found), 'out of', length(filter),
               filtername, 'were not found in the signature:',
               paste0(not_found, collapse = ", "))
  } else w <- NULL
  return(list(ids, w))
}
