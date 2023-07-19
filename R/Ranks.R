#' @title Ranks the signatures from most sensitive to least sensitive
#' @description  This function computes the beyondcell score's (BCS) statistics
#' of each signature and ranks them according to the residual's mean and the 
#' switch point. The signatures that pass the thresholds (see Details) are 
#' labelled as "TOP-LowSensitivity" or "TOP-HighSensitivity" (drugs to which 
#' all selected cells are least/most sensitive, respectively) or 
#' "TOP-Differential-LowSensitivity" and  "TOP-Differential-HighSensitivity" 
#' (drugs to which the cells of a specific cluster are differentially 
#' insensitive or sensitive when compared to the other clusters). These groups 
#' correspond to the quadrants drawn by \code{\link[beyondcell]{bc4Squares}}.
#' @name bcRanks
#' @import tidyverse
#' @importFrom data.table frank
#' @param bc \code{\link[beyondcell]{beyondcell}} object.
#' @param idents Name of the metadata column of interest. The signatures' ranks 
#' are computed for each level in \code{idents}.
#' @param extended If \code{extended = TRUE}, this function returns the switch
#' point, mean, median, standard deviation, variance, min, max, proportion of
#' \code{NaN} and residuals' mean per signature. If \code{extended = FALSE},
#' this function returns only the switch point, mean and residuals' mean.
#' @param resm.cutoff Numeric vector with the 2 desired residuals' mean cut-offs.
#' @param sp.cutoff Numeric vector with the 4 desired switch point cut-offs.
#' @details This function categorizes the signatures into 5 groups according to
#' the chosen residuals' mean (\code{resm.cutoff}) and switch point 
#' (\code{sp.cutoff}) cut-offs:
#' 
#' \itemize{
#' \item{TOP-LowSensitivity: signatures with a switch point > 
#' \code{max(sp.cutoff)} and a residuals' mean < \code{min(resm.cutoff)}.}
#' \item{TOP-HighSensitivity: signatures with a switch point < 
#' \code{min(sp.cutoff)} and a residuals' mean > \code{max(resm.cutoff)}.}
#' \item{TOP-Differential-LowSensitivity: signatures with a switch point 
#' between \code{sp.cutoff[2]} and \code{sp.cutoff[3]} and a residuals' mean <
#' \code{min(resm.cutoff)}.}
#' \item{TOP-Differential-HighSensitivity: signatures with a switch point 
#' between \code{sp.cutoff[2]} and \code{sp.cutoff[3]} and a residuals' mean >
#' \code{max(resm.cutoff)}.}
#' }
#' The signatures that do not meet the previous criteria are marked as
#' \code{NA}s.
#' Default residuals' mean cut-offs: first and last deciles; default switch 
#' point cut-offs: 0.1, 0.4, 0.6 and 0.9.
#' @return A \code{beyondcell} object with the results in a new entry of
#' \code{bc@@ranks}.
#' @examples
#' @export

bcRanks <- function(bc, idents = NULL, extended = TRUE, 
                    resm.cutoff = c(0.1, 0.9), 
                    sp.cutoff = c(0.1, 0.4, 0.6, 0.9)) {
  # Check that bc is a beyondcell object.
  if (class(bc) != "beyondcell") stop('bc must be a beyondcell object.')
  # Check idents.
  if (!is.null(idents)) {
    if (length(idents) != 1) {
      stop('Idents must be a single metadata column.')
    }
    if (idents %in% colnames(bc@meta.data)) {
      if (idents %in% names(bc@ranks)) {
        warning(paste0('$', idents, ' already exists in bc@ranks. ',
                       'Entry has been overwritten.'))
      }
      meta <- bc@meta.data[colnames(bc@normalized), idents, drop = TRUE]
    } else {
      stop('Idents not found.')
    }
  } else {
    stop("You must supply the name of a metadata column to group by.")
  }
  # Check extended.
  if (length(extended) != 1 | !is.logical(extended[1])) {
    stop('extended must be TRUE or FALSE.')
  }
  # Check resm.cutoff.
  if (length(resm.cutoff) != 2 | !is.numeric(resm.cutoff)) {
    stop('resm.cutoff must be a numeric vector of length 2.')
  }
  if (resm.cutoff[2] < resm.cutoff[1]) {
    warning(paste('Upper residuals\' mean cut-off is smaller than lower', 
                  'residuals\' mean cut-off. Sorting residuals\' mean cut-offs', 
                  'in increasing order.'))
    resm.cutoff <- sort(resm.cutoff, decreasing = FALSE)
  }
  # Check sp.cutoff.
  if (length(sp.cutoff) != 4 | !is.numeric(sp.cutoff)) {
    stop('sp.cutoff must be a numeric vector of length 4.')
  }
  if (any(sp.cutoff < 0 | sp.cutoff > 1)) {
    stop('sp.cutoff must contain 4 switch point values between 0 and 1.')
  }
  sorted.sp.cutoff <- sort(sp.cutoff, decreasing = FALSE)
  if (!identical(sp.cutoff, sorted.sp.cutoff)) {
    warning(paste('Sorting switch point cut-offs in increasing order.'))
    sp.cutoff <- sorted.sp.cutoff
  }
  # --- Code ---
  # Progress bar.
  pb <- txtProgressBar(min = 0, max = 100, style = 3, file = stderr())
  bins <- 10
  # Signatures in bc.
  sigs <- rownames(bc@normalized)
  # Cells in bc.
  cells <- colnames(bc@normalized)
  # Keep column to group by.
  meta <- bc@meta.data %>%
    tibble::rownames_to_column("cells") %>%
    dplyr::select(cells, all_of(idents)) %>%
    dplyr::rename(group.var := !!idents) %>%
    dplyr::mutate(group.var = factor(group.var)) %>%
    unique()
  lvls <- levels(meta$group.var)
  Sys.sleep(0.1)
  setTxtProgressBar(pb, value = 5)
  # Column to order by.
  order.col <- paste0("rank.", levels(meta$group.var)[1])
  # Final column order.
  if (extended) {
    cols.additional <- c("median", "sd", "variance", "min", "max", "prop.na")
  } else {
    cols.additional <- NULL
  }
  cols.stats <- c("rank", "switch.point", "mean", cols.additional, 
                  "residuals.mean", "group")
  cols.stats.level <- tidyr::expand_grid(lvls, cols.stats) %>%
    dplyr::mutate(col.name = paste(cols.stats, lvls, sep = ".")) %>%
    dplyr::pull(col.name)
  # Get switch points.
  sp <- data.frame(switch.point = bc@switch.point) %>%
    tibble::rownames_to_column("IDs")
  Sys.sleep(0.1)
  setTxtProgressBar(pb, value = 10)
  # Compute long normalized BCS.
  normalized.long <- bc@normalized %>%
    t() %>%
    as.data.frame() %>%
    tibble::rownames_to_column("cells") %>%
    tidyr::pivot_longer(cols = all_of(sigs), names_to = "IDs", 
                        values_to = "enrichment", values_drop_na = FALSE)
  # Add grouping information and switch point.
  normalized.long <- normalized.long %>%
    dplyr::inner_join(sp, by = "IDs") %>%
    dplyr::inner_join(meta, by = "cells")
  Sys.sleep(0.1)
  setTxtProgressBar(pb, value = 25)
  # Compute mean BCS 
  stats.long <- normalized.long %>%
    dplyr::group_by(IDs) %>%
    dplyr::mutate(mean = round(mean(enrichment, na.omit = TRUE), digits = 2)) %>%
    dplyr::ungroup()
  # Compute residual's mean per signature.
  invisible(capture.output(
    res.long <- normalized.long %>%
      dplyr::filter(!is.na(enrichment)) %>%
      dplyr::do(data.frame(., resid = residuals(lm(enrichment ~ group.var, data = .)))) %>%
      dplyr::group_by(IDs, group.var) %>%
      dplyr::mutate(residuals.mean = mean(resid, na.rm = TRUE)) %>%
      dplyr::ungroup() %>%
      dplyr::select(IDs, group.var, residuals.mean) %>%
      unique()
  ))
  stats.long <- stats.long %>%
    dplyr::left_join(res.long, by = c("IDs", "group.var"))
  Sys.sleep(0.1)
  setTxtProgressBar(pb, value = 45)
  # If extended == TRUE, compute the median, standard deviation, variance, min, 
  # max and proportion of NaNs per signature.
  if (extended) {
    stats.long <- stats.long %>%
      dplyr::group_by(IDs) %>%
      dplyr::mutate(median = round(median(enrichment, na.rm = TRUE), digits = 2),
                    sd = round(sd(enrichment, na.rm = TRUE), digits = 2),
                    variance = round(var(enrichment, na.rm = TRUE), digits = 2),
                    min = round(min(enrichment, na.rm = TRUE), digits = 2),
                    max = round(max(enrichment, na.rm = TRUE), digits = 2),
                    prop.na = round(sum(is.na(enrichment))/length(cells), 
                                    digits = 2)) %>%
      dplyr::ungroup()
  }
  Sys.sleep(0.1)
  setTxtProgressBar(pb, value = 50)
  # Residual's deciles.
  res.decil <- stats.long %>%
    dplyr::group_by(group.var) %>%
    dplyr::group_modify(~as.data.frame(t(quantile(.$residuals.mean, resm.cutoff, 
                                                  na.rm = TRUE)))) %>%
    dplyr::ungroup()
  colnames(res.decil)[2:3] <- c("Pmin", "Pmax")
  stats.long <- stats.long %>%
    dplyr::select(-cells, -enrichment) %>%
    unique() %>%
    dplyr::inner_join(res.decil, by = "group.var")
  Sys.sleep(0.1)
  setTxtProgressBar(pb, value = 75)
  # Group annotation.
  stats.long.annotated <- stats.long %>%
    dplyr::mutate(
      group = dplyr::case_when(switch.point < sp.cutoff[1] & 
                                 residuals.mean > Pmax ~ 
                                 "TOP-HighSensitivity",
                               switch.point > sp.cutoff[4] & 
                                 residuals.mean < Pmin ~ 
                                 "TOP-LowSensitivity",
                               switch.point > sp.cutoff[2] & 
                                 switch.point < sp.cutoff[3] & 
                                 residuals.mean < Pmin ~ 
                                 "TOP-Differential-LowSensitivity",
                               switch.point > sp.cutoff[2] &
                                 switch.point < sp.cutoff[3] & 
                                 residuals.mean > Pmax ~ 
                                 "TOP-Differential-HighSensitivity",
                               TRUE ~ NA_character_))
  Sys.sleep(0.1)
  setTxtProgressBar(pb, value = 80)
  # Order.
  rank <- stats.long.annotated %>%
    dplyr::mutate(in.range = switch.point > sp.cutoff[2] & 
                    switch.point < sp.cutoff[3],
                  sp.rank = switch.point * as.numeric(in.range)) %>%
    dplyr::select(IDs, group.var, sp.rank, residuals.mean, in.range) %>%
    unique() %>%
    dplyr::group_split(group.var)
  rank <- lapply(rank, FUN = function(x) {
    dt <- data.table::as.data.table(x)
    dt[, rank := data.table::frank(dt, -sp.rank, -residuals.mean, 
                                   ties.method = "dense")]
    return(dt)
  }) %>%
    dplyr::bind_rows() %>%
    dplyr::mutate(rank = dplyr::if_else(in.range, rank, NA_integer_)) %>%
    dplyr::select(IDs, group.var, rank) %>%
    unique()
  stats.long.ranked <- stats.long.annotated %>%
    dplyr::inner_join(rank, by = c("IDs", "group.var"))
  Sys.sleep(0.1)
  setTxtProgressBar(pb, value = 85)
  # Pivot wider
  final.stats <- stats.long.ranked %>%
    dplyr::select(IDs, group.var, all_of(cols.stats)) %>%
    unique() %>%
    tidyr::pivot_wider(names_from = group.var, values_from = all_of(cols.stats),
                       names_sep = ".")
  Sys.sleep(0.1)
  setTxtProgressBar(pb, value = 90)
  # Add Drug name and MoA to final.stats.
  info <- drugInfo$IDs %>%
    dplyr::filter(IDs %in% final.stats$IDs) %>%
    dplyr::select(IDs, preferred.drug.names, studies) %>%
    dplyr::left_join(y = drugInfo$MoAs[, c("IDs", "MoAs")], by = "IDs",
                     relationship = "many-to-many") %>%
    dplyr::left_join(y = drugInfo$Targets, by = "IDs",
                     relationship = "many-to-many") %>%
    dplyr::left_join(y = drugInfo$Synonyms, by = "IDs",
                     relationship = "many-to-many")
  if (dim(info)[1] > 0) {
    info <- aggregate(.~ IDs, data = info, na.action = NULL, FUN = function(x) {
      paste(na.omit(unique(x)), collapse = "; ")
    })
    cols.druginfo <- c("drugs", "preferred.drug.names", "MoAs", "targets", 
                       "studies")
  } else {
    info <- data.frame(IDs = rownames(bc@normalized))
    cols.druginfo <- NULL
  }
  final.stats <- final.stats %>%
    dplyr::inner_join(info, by = "IDs") %>%
    tibble::column_to_rownames("IDs") %>%
    unique()
  Sys.sleep(0.1)
  setTxtProgressBar(pb, value = 95)
  # Order by rank and reorder columns.
  final.stats <- final.stats[order(final.stats[, order.col], decreasing = FALSE),
                             c(cols.druginfo, cols.stats.level)]
  # Add to beyondcell object.
  bc@ranks[[idents]] <- final.stats
  Sys.sleep(0.1)
  setTxtProgressBar(pb, value = 100)
  return(bc)
}

#' @title Returns the first/last n ranked signatures
#' @description  This function returns the top/bottom \code{n} signatures ranked
#' by \code{\link[beyondcell]{bcRanks}}. If the rank has not been previously
#' computed, \code{rankSigs} performs the ranking itself.
#' @name rankSigs
#' @param bc \code{\link[beyondcell]{beyondcell}} object.
#' @param idents Name of the metadata column of interest. If
#' \code{idents = NULL}, the function uses the general rank computed with all
#' cells.
#' @param cond Level of \code{idents} to rank by the output vector. If
#' \code{idents = NULL}, this parameter is deprecated.
#' @param n Number of signatures to return in the output vector.
#' @param decreasing Logical. Return the top \code{n} signatures (default) or
#' the bottom \code{n} signatures (\code{decreasing = FALSE})?.
#' @return An ordered vector with the signature's names.
#' @examples
#' @export

rankSigs <- function(bc, idents = NULL, cond = NULL, n = 10,
                     decreasing = TRUE) {
  # --- Checks ---
  # Check that bc is a beyondcell object.
  if (class(bc) != "beyondcell") stop('bc must be a beyondcell object.')
  # Check idents.
  if (!is.null(idents)) {
    if (length(idents) != 1) {
      stop('Idents must be a single metadata column.')
    }
    if (!idents %in% colnames(bc@meta.data)) {
      stop('Idents not found.')
    }
    if (is.null(cond[1])) {
      stop('Invalid cond.')
    }
    meta <- idents
  } else {
    meta <- "general"
    if (!is.null(cond[1])) {
      warning('idents not specified, cond is deprecated.')
      cond <- NULL
    }
  }
  # Check cond.
  if (!is.null(cond)) {
    if (length(cond) != 1) {
      stop('cond must be a single idents level')
    }
    if (!cond %in% levels(as.factor(bc@meta.data[, idents]))) {
      stop(paste0(cond, ' is not a level of ', idents, '.'))
    }
  }
  # Check n.
  if (length(n) != 1 | (!is.numeric(n) & !is.character(n))) {
    stop('n must be a single number or "all".')
  }
  if (is.numeric(n) & (n < 1 | n%%1 != 0)) stop('n must be an integer > 0.')
  else if (is.character(n)) {
    if (n == "all") n <- nrow(bc@normalized)
    else stop('To select all signatures, please set n = "all".')
  }
  # Check decreasing.
  if (length(decreasing) != 1 | !is.logical(decreasing)) {
    stop('decreasimg must be TRUE or FALSE.')
  }
  # --- Code ---
  # If ranks have not been computed, compute them now.
  if (!meta %in% names(bc@ranks)) {
    message('Computing ranks...')
    bc <- bcRanks(bc, idents = idents, extended = FALSE)
  }
  # Get ranks for the specified idents.
  df <- bc@ranks[[meta]]
  # Ranks to select.
  if (decreasing == TRUE) {
    idx <- 1:n
  } else {
    idx <- nrow(bc@normalized):(nrow(bc@normalized) - n + 1)
  }
  # Return signatures whose rank == idx.
  order.col <- ifelse(test = meta == "general", yes = "rank",
                      no = paste0("rank.", cond))
  ordered.df <- df[order(df[, order.col], decreasing = FALSE), ]
  sigs <- ordered.df[idx, "Name"]
  return(sigs)
}
