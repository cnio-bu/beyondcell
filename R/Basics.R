#' @title Substracts to the first element of a vector the rest of elements
#' @description This function substracts to the first element of a numeric
#' vector (\code{x[1]}) the rest of elements of the same vector.
#' (\code{x[2:length(x)]}).
#' @name minus
#' @param x Numeric vector.
#' @param na.rm (From \code{base}) Logical. Should missing values (including
#' \code{NaN}) be removed?
#' @return The result of the substraction.
#' @examples
#' @export

minus <- function(x, na.rm = FALSE) {
  # --- Checks ---
  # Check x.
  if (!is.numeric(x)) stop('x must be numeric.')
  # Check na.rm.
  if (length(na.rm) != 1 | !is.logical(na.rm)) {
    stop('na.rm must be TRUE or FALSE.')
  }
  # --- Code ---
  # If x is a single number, append a zero so we can run the next step.
  if (length(x) == 1) x <- c(x, 0)
  # Substract to the first element of x the rest of elements.
  out <- sum(x[1], na.rm = na.rm) - sum(x[2:length(x)], na.rm = na.rm)
  return(out)
}

#' @title Computes the column substraction
#' @description This function substracts to the first element of each column of
#' a rectangular object (\code{x[1, n]}) the rest of elements of the same column
#' (\code{x[2:length(x), n]}).
#' @name colMinus
#' @param x A matrix or a dataframe.
#' @param na.rm (From \code{base}) Logical. Should missing values (including
#' \code{NaN}) from rows \code{2:length(x)} be omitted from the calculations?
#' @return A numeric rectangular object with the result of the substraction.
#' @examples
#' @export

colMinus <- function(x, na.rm = FALSE) {
  # --- Checks ---
  # Check x.
  if (!is.matrix(x) & !is.data.frame(x)) {
    stop('x must be a matrix or a data.frame.')
  }
  # Check na.rm.
  if (length(na.rm) != 1 | !is.logical(na.rm)) {
    stop('na.rm must be TRUE or FALSE.')
  }
  # --- Code ---
  # If x has a single row, append a row of zeros so we can run the next step.
  if (dim(x)[1] == 1) x <- rbind(x, rep(0, times = ncol(x)))
  # Substract to the first row of x the rest of rows.
  first.row <- x[1, , drop = FALSE]
  out <- first.row - colSums(x[2:nrow(x), , drop = FALSE], na.rm = na.rm)
  return(out)
}

#' @title Computes the mean, the median and the sd of a vector
#' @description This function computes the mean, the median and the sd of a
#' vector.
#' @name Mean.Med.SD
#' @param x Numeric vector.
#' @param na.rm (From base) Logical. Should missing values (including NaN) be
#' removed?
#' @return A named numeric vector with the mean, median and sd of \code{x}.
#' @examples
#' @export

Mean.Med.SD <- function(x) {
  # --- Checks ---
  # Check that x is a numeric vector.
  if (length(x) <= 1 | !is.numeric(x)) {
    stop('x must be a numeric vector.')
  }
  # --- Code ---
  stats.mean <- mean(x, na.rm = TRUE)
  stats.median <- median(x, na.rm = TRUE)
  stats.sd <- sd(x, na.rm = TRUE)
  return(setNames(c(stats.mean, stats.median, stats.sd),
                  c("mean", "median", "sd")))
}

#' @title Returns the IDs that match the specified filter's values
#' @description This function subsets \code{df} to select only the entries that
#' match the specified \code{filter}'s \code{values} and returns the
#' corresponding IDs.
#' @name GetIDs
#' @param values User-supplied filtering vector.
#' @param filter Column name or number to subset by.
#' @return A vector with the IDs that match the \code{filter}'s values.
#' @export

GetIDs <- function(values, filter) {
  # --- Checks ---
  # Check values.
  if (length(values) < 1 | !is.character(values)) {
    stop('values must be a character vector.')
  }
  # Check filter.
  if (length(filter) != 1) {
    stop('You must specify a single filter.')
  }
  all_filters <- c("Synonyms", "IDs", "MoAs", "Targets", "IDs")
  names(all_filters) <- c("drugs", "IDs", "MoAs", "targets", "studies")

  df <- drugInfo[[all_filters[filter]]]
  
  # --- Code ---
  upper.values <- toupper(values)
  selected <- subset(df, subset = toupper(df[[filter]]) %in% upper.values)
  if (filter == "drugs" & ".drug.names" %in% colnames(df)) {
    synonyms <- subset(df, subset = toupper(df[["preferred.drug.names"]]) %in%
                         unique(toupper(selected[["preferred.drug.names"]])))
    selected <- unique(rbind(selected, synonyms))
  }
  ids <- unique(selected$IDs)
  not.found <- values[!(upper.values %in% toupper(df[[filter]]))]
  if (all(values %in% not.found)) {
    stop('No sig ID was found for any of the elements in values.')
  } else if (length(not.found) > 0) {
    filtername <- gsub(pattern = '"', replacement = '',
                       x = deparse(substitute(filter)))
    warning(paste0('sig IDs were not found for ', length(not.found), ' out of ',
                   length(values), " ", filtername, ': ',
                   paste0(not.found, collapse = ", ")))
  }
  return(ids)
}

#' @title Returns a dataframe with information about the input drugs
#' @description This function searches the input drugs in the pre-loaded
#' \code{beyondcell} matrices and returns a dataframe with drug information,
#' including drug synonyms and MoAs.
#' @name FindDrugs
#' @importFrom dplyr left_join
#' @param bc \code{\link[beyondcell]{beyondcell}} object.
#' @param x A character vector with drug names and/or sig IDs.
#' @param na.rm Logical. Should \code{x} entries with no available drug 
#' information be removed from the final output?
#' @details The output \code{data.frame} has the following columns:
#' \itemize{
#' \item{\code{original.names}}: Input drug names.
#' \item{\code{bc.names}}: Drug names used in \code{bc}.
#' \item{\code{preferred.drug.names}}: Standard drug names.
#' \item{\code{drugs}}: Other drug names.
#' \item{\code{IDs}}: Signature IDs.
#' \item{\code{preferred.and.sigs}}: \code{preferred.drug.names} (or alternatively
#' \code{bc.names}) and \code{IDs}. Used as title in \code{beyondcell} plots.
#' \item{\code{MoAs}}: Mechanism(s) of action.
#' }
#' @return A \code{data.frame}.
#' @examples
#' @export

FindDrugs <- function(bc, x, na.rm = TRUE) {
  # --- Checks ---
  # Check that bc is a beyondcell object.
  if (class(bc) != "beyondcell") stop('bc must be a beyondcell object.')
  # Check x.
  if (!is.character(x)) stop('x must be a character vector.')
  # Check na.rm.
  if (length(na.rm) != 1 | !is.logical(na.rm)) {
    stop('na.rm must be TRUE or FALSE.')
  }
  # --- Code ---
  # workaround to expose a facade of the old database to this function
  # instead of the list of data.frames that we have now.
  old_db <- drugInfo$IDs %>%
    dplyr::select(IDs, preferred.drug.names, studies) %>%
    dplyr::left_join(y = drugInfo$MoAs[, c("IDs", "MoAs")], by = "IDs") %>%
    dplyr::left_join(y = drugInfo$Targets, by = "IDs") %>%
    dplyr::left_join(y = drugInfo$Synonyms, by = "IDs") %>%
    dplyr::mutate(
      sources = studies
    ) %>%
    as.data.frame()
  
  # bc signatures.
  sigs <- rownames(bc@normalized)
  # Match x with bc signatures and get the indexes of matching elements.
  indexes <- lapply(x, function(y) {
    idx <- match(toupper(y), table = toupper(sigs), nomatch = 0)
    if (idx == 0) {
      idx <- unique(match(old_db$IDs[old_db$drugs == toupper(y)],
                          table = sigs))
    }
    return(idx[!is.na(idx)])
  })
  # Original names (x) and bc names (sigs).
  df <- data.frame(original.names = unlist(sapply(seq_along(x), function(i) {
    rep(x[i], times = length(indexes[[i]]))
  })), IDs = unlist(sapply(indexes, function(z) sigs[z])))
  df.not.found <- !(x %in% df$original.names)
  if (any(df.not.found)) {
    empty.df <- data.frame(original.names = x[df.not.found],
                           IDs = rep(NA, sum(df.not.found)))
    df <- rbind(df, empty.df)
  }
  # Get the names and pathways of the selected signatures.
  info <- subset(old_db, subset = IDs %in% df$IDs)
  if (all(dim(info) != 0)) {
    info <- aggregate(.~ IDs, data = info, na.action = NULL, FUN = function(w) {
      paste(na.omit(unique(w)), collapse = ", ")
    })
  }
  info.not.found <- !(df$IDs %in% old_db$IDs)
  if (any(info.not.found)) {
    empty.info <- matrix(rep(NA,
                             times = sum(info.not.found)*ncol(old_db)),
                         ncol = ncol(old_db),
                         dimnames = list(1:sum(info.not.found),
                                         colnames(old_db)))
    info <- rbind(info, as.data.frame(empty.info))
  }
  # Merge df and info.
  df <- unique(merge(df, info[, c("IDs", "drugs", "preferred.drug.names",
                                  "MoAs")], by = "IDs", all.x = TRUE))
  # Add bc.names column and remove names that are not sig IDs from sig_id
  # column.
  df$bc.names <- df$IDs
  df$IDs[!startsWith(df$IDs, prefix = "sig-")] <- NA
  # Create preferred.and.sigs column: Preferred name and sig_id.
  df$preferred.and.sigs <- sapply(1:nrow(df), function(j) {
    return(ifelse(test = !is.na(df$preferred.drug.names[j]),
                  yes = paste0(df$preferred.drug.names[j],
                               paste0(" (", df$IDs[j], ")")),
                  no = df$bc.names[j]))
  })
  # Reorder df.
  rows <- unlist(lapply(x, function(entry) which(df$original.names == entry)))
  cols <- c("original.names", "bc.names", "preferred.drug.names", "drugs", "IDs",
            "preferred.and.sigs", "MoAs")
  df <- df[rows, cols]
  # If na.rm = TRUE, remove rows with NAs in "preferred.drug.names" and "drugs" 
  # fields.
  if (na.rm) df <- df[rowSums(is.na(df[, 3:4])) < 2, ]
  return(df)
}
