#' @title Ranks the signatures from most sensitive to least sensitive
#' @description  This function computes the beyondcell score's (BCS) statistics
#' of each signature and ranks them according to the switch point and 
#' residual's mean.
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
#' @return A \code{beyondcell} object with the results in a new entry of
#' \code{bc@@ranks}.
#' @examples
#' @export

bcRanks <- function(bc, idents = NULL, extended = TRUE) {
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
  # --- Code ---
  # Signatures in bc.
  sigs <- rownames(bc@normalized)
  # Cells in bc.
  cells <- colnames(bc@normalized)
  # Keep column to group by.
  meta <- bc@meta.data %>%
    rownames_to_column("cells") %>%
    select(cells, all_of(idents)) %>%
    rename(group.var := !!idents) %>%
    mutate(group.var = factor(group.var)) %>%
    unique()
  lvls <- levels(meta$group.var)
  # Column to order by.
  order.col <- paste0("rank.", levels(meta$group.var)[1])
  # Final column order.
  if (extended) {
    cols.additional <- c("median", "sd", "variance", "min", "max", "prop.na")
  } else {
    cols.additional <- NULL
  }
  cols.druginfo <- c("drugs", "preferred.drug.names", "MoAs", "targets", 
                     "studies")
  cols.stats <- c("rank", "switch.point", "mean", cols.additional, 
                  "residuals.mean", "group")
  cols.stats.level <- expand_grid(lvls, cols.stats) %>%
    mutate(col.name = paste(cols.stats, lvls, sep = ".")) %>%
    pull(col.name)
  # Get switch points.
  sp <- data.frame(switch.point = bc@switch.point) %>%
    rownames_to_column("IDs")
  # Compute long normalized BCS.
  normalized.long <- bc@normalized %>%
    t() %>%
    as.data.frame() %>%
    rownames_to_column("cells") %>%
    pivot_longer(cols = all_of(sigs), names_to = "IDs", 
                 values_to = "enrichment", values_drop_na = FALSE)
  # Add grouping information and switch point.
  normalized.long <- normalized.long %>%
    dplyr::inner_join(sp, by = "IDs") %>%
    dplyr::inner_join(meta, by = "cells")
  # Compute mean BCS and residual's mean per signature.
  stats.long <- normalized.long %>%
    group_by(IDs) %>%
    mutate(mean = round(mean(enrichment), 2)) %>%
    do(data.frame(., resid = residuals(lm(enrichment ~ group.var, data = .)))) %>%
    group_by(IDs, group.var) %>%
    mutate(residuals.mean = mean(resid, na.rm = TRUE)) %>%
    ungroup()
  # If extended == TRUE, compute the median, standard deviation, variance, min, 
  # max and proportion of NaNs per signature.
  if (extended) {
    stats.long <- stats.long %>%
      group_by(IDs) %>%
      mutate(median = round(median(enrichment, na.rm = TRUE), digits = 2),
             sd = round(sd(enrichment, na.rm = TRUE), digits = 2),
             variance = round(var(enrichment, na.rm = TRUE), digits = 2),
             min = round(min(enrichment, na.rm = TRUE), digits = 2),
             max = round(max(enrichment, na.rm = TRUE), digits = 2),
             prop.na = round(sum(is.na(enrichment))/length(cells), digits = 2)) %>%
      ungroup()
  }
  # Residual's deciles.
  res.decil <- stats.long %>%
    group_by(group.var) %>%
    group_modify(~as.data.frame(t(quantile(.$residuals.mean, c(0.1, 0.9))))) %>%
    ungroup()
  colnames(res.decil)[2:3] <- c("P10", "P90")
  stats.long <- stats.long %>%
    inner_join(res.decil, by = "group.var") %>%
    select(-cells, -enrichment, -resid) %>%
    unique()
  # Group annotation.
  stats.long.annotated <- stats.long %>%
    mutate(group = case_when(switch.point < 0.1 & residuals.mean > P90 ~ 
                               "TOP-HighSensitivity",
                             switch.point > 0.9 & residuals.mean < P10 ~ 
                               "TOP-LowSensitivity",
                             switch.point > 0.4 & switch.point < 0.6 &
                               residuals.mean < P10 ~ 
                               "TOP-Differential-LowSensitivity",
                             switch.point > 0.4 & switch.point < 0.6 & 
                               residuals.mean > P90 ~ 
                               "TOP-Differential-HighSensitivity",
                             TRUE ~ NA_character_))
  # Order.
  rank <- stats.long.annotated %>%
    mutate(in.range = switch.point > 0.4 & switch.point < 0.6,
           sp.rank = switch.point * as.numeric(in.range)) %>%
    select(IDs, group.var, sp.rank, residuals.mean, in.range) %>%
    unique() %>%
    group_split(group.var)
  rank <- lapply(rank, FUN = function(x) {
    dt <- as.data.table(x)
    dt[, rank := data.table::frank(dt, -sp.rank, -residuals.mean, 
                                   ties.method = "dense")]
    return(dt)
  }) %>%
    bind_rows() %>%
    mutate(rank = if_else(in.range, rank, NA_integer_)) %>%
    select(IDs, group.var, rank) %>%
    unique()
  stats.long.ranked <- stats.long.annotated %>%
    inner_join(rank, by = c("IDs", "group.var"))
  # Pivot wider
  final.stats <- stats.long.ranked %>%
    select(IDs, group.var, all_of(cols.stats)) %>%
    unique() %>%
    pivot_wider(names_from = group.var, values_from = all_of(cols.stats),
                names_sep = ".")
  # Add Drug name and MoA to final.stats.
  info <- drugInfo$IDs %>%
    filter(IDs %in% final.stats$IDs) %>%
    select(IDs, preferred.drug.names, studies) %>%
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
  }
  final.stats <- final.stats %>%
    inner_join(info, by = "IDs") %>%
    column_to_rownames("IDs") %>%
    unique()
  # Order by rank and reorder columns.
  final.stats <- final.stats[order(final.stats[, order.col], decreasing = FALSE),
                             c(cols.druginfo, cols.stats.level)]
  # Add to beyondcell object.
  bc@ranks[[idents]] <- final.stats
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
