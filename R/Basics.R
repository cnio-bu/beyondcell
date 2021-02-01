#' @title Substracts to the first element of a vector the rest of elements
#' @description This function substracts to the first element of a numerical
#' vector (\code{x[1]}) the rest of elements of the same vector.
#' (\code{x[2:length(x)]}).
#' @name minus
#' @param x Numerical vector
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

#' @title Substracts to the first row of a rectangular object the rest of rows
#' @description This function substracts to the first row of a numerical
#' rectangular object (\code{x[1, ]}) the rest of rows of the same rectangular
#' object (\code{x[2:length(x), ]}).
#' @name colMinus
#' @param x A matrix or a data frame.
#' @param na.rm (From \code{base}) Logical. Should missing values (including
#' \code{NaN}) from rows 2:length(x) be omitted from the calculations?
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
#' vector removing \code{NA} values.
#' @name Mean.Med.SD
#' @param x A numeric vector for which to compute the statistics.
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

#' @title Computes bcRanks' statistics and ranks
#' @description This function computes the bcscores' statistics and ranks
#' returned by \code{\link[bcRnaks]{bcRanks}}.
#' @name GetStatistics
#' @param bc \code{\link[beyondcell]{beyondcell}} object.
#' @param signatures Vector with the names of the signatures of interest.
#' @param cells Vector with the names of the cells of interest.
#' @param pb Progress bar.
#' @param total Number of iterations to complete the progress bar.
#' @param i Iteration number; used for increasing the progress bar.
#' @param n.rows Number of signatures; used for increasing the progress bar.
#' @param extended If \code{extended = TRUE}, this function returns the switch
#' point, mean, median, sd, variance, min, max, proportion of \code{NaN} and
#' residuals' mean per signature. If \code{extended = FALSE}, this function
#' returns only the switch point, mean and residuals' mean.
#'
#' @return A \code{data.frame} with the statistics and ranks per signature.
#' @examples
#' @export

GetStatistics <- function(bc, signatures, cells, pb, total, i, n.rows,
                          extended) {
  # --- Checks ---
  # Check that bc is a beyondcell object.
  if (class(bc) != "beyondcell") stop('bc must be a beyondcell object.')
  # Check signatures.
  in.signatures <- !is.null(signatures) & signatures %in% rownames(bc@normalized)
  if (all(!in.signatures)) {
    stop('None of the specified signatures were found.')
  } else if (any(!in.signatures)) {
    warning(paste0('These signatures were not found in bc: ',
                   paste0(signatures[!in.signatures], collapse = ", "), '.'))
  }
  # Check cells.
  in.cells <- !is.null(cells) & cells %in% colnames(bc@normalized)
  if (all(!in.cells)) {
    stop('None of the specified cells were found.')
  } else if (any(!in.cells)) {
    warning(paste0('These cells were not found in bc: ',
                   paste0(cells[!in.cells], collapse = ", "), '.'))
  }
  # Check pb.
  if (class(pb) != "txtProgressBar") stop('pb must be a txtProgressBar object.')
  # Check total.
  if (length(total) != 1 | total[1] <= 0 | total[1]%%1 != 0) {
    stop('total must be a positive integer.')
  }
  # Check i.
  if (length(i) != 1 | i[1] <= 0 | i[1]%%1 != 0) {
    stop('i must be a positive integer.')
  }
  # Check n.rows.
  if (length(n.rows) != 1 | n.rows[1] <= 0 | n.rows[1]%%1 != 0) {
    stop('n.rows must be a positive integer.')
  }
  # Check extended.
  if (length(extended) != 1 | !is.logical(extended[1])) {
    stop('extended must be TRUE or FALSE.')
  }
  # --- Code ---
  # Signatures and cells.
  signatures <- signatures[in.signatures]
  cells <- cells[in.cells]
  # Bins: Number of iterations to complete each 1% of the progress bar.
  bins <- ceiling(total / 100)
  # Switch points per signature.
  switch.p <- bc@switch.point[signatures]
  if (extended) {
    # Dataframe with mean, median and sd per signature.
    data <- cbind(seq_len(n.rows) + (n.rows * 6) * (i - 1),
                  bc@data[signatures, cells])
    mean.med.sd <- as.data.frame(t(apply(data, 1, function(u) {
      mms <- round(Mean.Med.SD(u[-1]), digits = 2)
      ### Update the progress bar.
      if (u[1]%%bins == 0) {
        Sys.sleep(0.1)
        setTxtProgressBar(pb, value = u[1])
      }
      return(mms)
    })))
    # Variance per signature.
    data[, 1] <- data[, 1] + n.rows
    variance.bcscore <- apply(data, 1, function(v) {
      variance <- var(v[-1], na.rm = TRUE)
      ### Update the progress bar.
      if (v[1]%%bins == 0) {
        Sys.sleep(0.1)
        setTxtProgressBar(pb, value = v[1])
      }
      return(variance)
    })
    # Min normalized BCS per signature.
    data[, 1] <- data[, 1] + n.rows
    min.bcscore <- apply(data, 1, function(w) {
      min.bcs <- min(w[-1], na.rm = TRUE)
      ### Update the progress bar.
      if (w[1]%%bins == 0) {
        Sys.sleep(0.1)
        setTxtProgressBar(pb, value = w[1])
      }
      return(min.bcs)
    })
    # Max normalized BCS per signature.
    data[, 1] <- data[, 1] + n.rows
    max.bcscore <- apply(data, 1, function(x) {
      max.bcs <- max(x[-1], na.rm = TRUE)
      ### Update the progress bar.
      if (x[1]%%bins == 0) {
        Sys.sleep(0.1)
        setTxtProgressBar(pb, value = x[1])
      }
      return(max.bcs)
    })
    # NA proportion per signature.
    data[, 1] <- data[, 1] + n.rows
    prop.na <- apply(data, 1, function(y) {
      nas <- round(sum(is.na(y[-1]))/length(y[-1]), digits = 2)
      ### Update the progress bar.
      if (y[1]%%bins == 0) {
        Sys.sleep(0.1)
        setTxtProgressBar(pb, value = y[1])
      }
      return(nas)
    })
    # Residuals.
    normalized <- cbind(seq_len(n.rows) + (n.rows * i * 5) + (n.rows * (i - 1)),
                        bc@normalized[signatures, cells])
    # Create dataframe.
    stats <- data.frame(switch.point = switch.p, mean = mean.med.sd$mean,
                        median = mean.med.sd$median, sd = mean.med.sd$sd,
                        variance = variance.bcscore, min = min.bcscore,
                        max = max.bcscore, prop.na = prop.na,
                        row.names = signatures)
  } else {
    # Mean BCS per signature.
    mean.bc <- round(rowMeans(bc@data[signatures, cells], na.rm = TRUE),
                     digits = 2)
    # Residuals.
    normalized <- cbind(seq_len(n.rows) + (n.rows * (i - 1)),
                        bc@normalized[signatures, cells])
    # Create dataframe.
    stats <- data.frame(switch.point = switch.p, mean = mean.bc,
                        row.names = signatures)
  }
  # Regression residuals.
  resid <- apply(normalized, 1, function(z) {
    res <- round(mean(z[-1], na.rm = TRUE), digits = 2)
    ### Update the progress bar.
    if (z[1]%%bins == 0 | z[1] == total) {
      Sys.sleep(0.1)
      setTxtProgressBar(pb, value = z[1])
    }
    return(res)
  })
  # Update dataframe.
  stats <- cbind(stats, data.frame(residuals = resid, row.names = signatures))
  return(stats)
}
