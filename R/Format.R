#' @title Returns a vector of colors to use in \code{center_scale_colour_stepsn}
#' @description This function gets a color palette and returns a vector with 5
#' color values to be used in
#' \code{\link[center_scale_colour_stepsn]{center_scale_colour_stepsn}}.
#' @name get_colour_stepsn
#' @importFrom viridis viridis
#' @importrom RColorBrewer brewer.pal
#' @param colorscale Either a \code{viridis}, \code{RColorBrewer} or a custom
#' palette of 3 colors (low, medium and high). If \code{colorscale = NULL}
#' (default), the function returns \code{beyondcell}'s own palette.
#' @return A vector with 5 color values.
#' @examples
#' @export

get_colour_stepsn <- function(colorscale = NULL) {
  # --- Checks and Code ---
  default <- c("#1D61F2", "#98B9FF", "#F7F7F7", "#FF9CBB", "#DA0078")
  # Check colorscale and get its value.
  if (is.null(colorscale)) {
    colors <- default
  } else {
    ### Try catchs:
    ### Is colorscale a viridis palette?
    guess.colors <- tryCatch(viridis::viridis(11, option = colorscale),
                             error = function(cond) cond,
                             warning = function(cond) cond)
    if (inherits(guess.colors, "error") | inherits(guess.colors, "warning")) {
      ### Is colorscale an RColorBrewer palette?
      guess.colors <- tryCatch(suppressWarnings(
        RColorBrewer::brewer.pal(12, name = colorscale)),
        error = function(cond) cond)
      if (inherits(guess.colors, "error")) {
        ### Is colorscale any other palette?
        guess.colors <- tryCatch(scale_colour_stepsn(colours = colorscale),
                                 error = function(cond) cond)
        if (inherits(guess.colors, "error")) {
          ### If not, set default colorscale.
          warning('Colorscale not found. Settig default colorscale.')
          colors <- default
        } else {
          ### If colorscale contains less than 3 values, set default colorscale.
          len.colors <- length(colorscale)
          if (len.colors < 3) {
            warning(paste('Colorscale too short. It must contain 3 colours:',
                          'high, medium and low. Settig default colorscale...'))
            colors <- default
          } else {
            ### Else, construct a scale with 5 colours: the first, middle and
            ### last values in colorscale and 2 intermediate colours between
            ### them (color.low and color.high, computed with colorRampPalette).
            color.middle <- colorscale[ceiling(len.colors/2)]
            color.low <- colorRampPalette(colors = c(colorscale[1], color.middle),
                                          space = "Lab")(3)[2]
            color.high <- colorRampPalette(colors = c(colorscale[3], color.middle),
                                           space = "Lab")(3)[2]
            colors <- c(colorscale[1], color.low, color.middle, color.high,
                        colorscale[len.colors])
            if (len.colors > 3) {
              warning(paste('Colorscale too long. It must contain 3 colours:',
                            'high, medium and low. Colours chosen:',
                            paste0(colors[c(1, 3, 5)], collapse = ", ")))
            }
          }
        }
        ### If colorscale is an RColorBrewer palette, subset 5 values to create
        ### the final palette.
      } else {
        len.guess <- length(guess.colors)
        idx.middle <- ceiling(len.guess/2)
        colors <- guess.colors[c(1, idx.middle - 1, idx.middle,
                                 idx.middle + 1, len.guess)]
      }
      ### If colorscale is a viridis palette, subset 5 values to create the final
      ### palette.
    } else {
      colors <- guess.colors[c(1, 5, 6, 7, 11)]
    }
  }
  return(colors)
}

#' @title Creates a centered diverging binned colour gradient
#' @description This function creates a diverging binned colour gradient
#' (low-mid-high) centered around \code{center}.
#' @name center_scale_colour_stepsn
#' @import ggplot2
#' @import scales
#' @param x A numeric vector. Can contain \code{NA}s.
#' @param colorscale A vector with 5 colors which can be obtained using
#' \code{\link[get_colour_stepsn]{get_colour_stepsn}}.
#' @param alpha Transparency level between 0 (not transparent) and 1 (fully
#' transparent).
#' @param na.value Color to use for missing values.
#' @param limits Vector with the desired limits.
#' @param center A single number indicating the center of the \code{colorscale}.
#' If \code{center = NULL} (default), the center is set to the middle point of
#' \code{x}.
#' @param breaks A single number indicating the break size of the
#' \code{colorscale}. Alternatively, it can be a vector with the desired breaks
#' (which don't have to be symmetric or equally distributed). If \code{center}
#' is a vector of two numbers, \code{breaks} are computed using the difference
#' between them.
#' @return A centered diverging binned colour gradient that can be use to color
#' \code{\link[ggplot2]{ggplot2}} objects.
#' @examples
#' @export

center_scale_colour_stepsn <- function(x, colorscale, alpha = 0.7,
                                       na.value = "grey50", limits = c(NA, NA),
                                       center = NULL, breaks = 0.1) {
  # --- Checks ---
  # Check x.
  if (!is.numeric(x)) {
    stop('x must be a numeric vector.')
  }
  range.values <- pretty(x)
  # Check colorscale.
  if (length(colorscale) != 5 |
      !tryCatch(is.matrix(col2rgb(colorscale)), error = function(e) FALSE)) {
    stop('colorscale must contain exactly 5 colours.')
  }
  # Check alpha.
  if (length(alpha) != 1 | alpha[1] < 0 | alpha[1] > 1) {
    stop('alpha must be a positive number between 0 and 1.')
  }
  # Check na.value.
  if (!tryCatch(is.matrix(col2rgb(na.value)), error = function(e) FALSE)) {
    stop('na.value is not a colour.')
  }
  # Check limits.
  if (length(limits) != 2) {
    stop('limits must be a vector of length 2.')
  }
  na.limits <- is.na(limits)
  if (length(limits[!na.limits]) > 0 & !is.numeric(limits[!na.limits])) {
    stop('limits must be numeric or NAs.')
  }
  # If some limits are NAs, compute them.
  if (any(na.limits)) {
    limits[na.limits] <- c(min(range.values), max(range.values))[na.limits]
  }
  # If limits are not sorted, sort them.
  if (limits[2] < limits[1]) {
    warning(paste('Upper limit is smaller than lower limit.',
                  'Sorting limits in increasing order.'))
    limits <- sort(limits, decreasing = FALSE)
  }
  # Check center.
  if (!is.null(center)) {
    if (length(center)!= 1| !is.numeric(center)) {
      stop('center must be a single number.')
    }
    if (center < limits[1] | center > limits[2]) {
      stop(paste('center =', center, 'outside of limits =',
                 paste0(limits, collapse = ", ")))
    }
    # If center = NULL, set center to middle point in range.values.
  } else {
    len.range <- length(range.values)
    ### If len.range is odd, get the middle point.
    if (len.range%%2 == 1) {
      center <- range.values[ceiling(len.range/2)]
      ### If len.range is even, get the two middle points and do the mean.
    } else if (len.range%%2 == 0) {
      center <- round(sum(range.values[(len.range/2):((len.range/2)+1)])/2,
                      digits = 2)
    }
  }
  # Check breaks.
  if (!is.numeric(breaks)) {
    stop('breaks must be numeric.')
  }
  # If breaks is a single number...
  if (length(breaks) == 1) {
    if (breaks > abs(limits[1] - limits[2])) {
      stop('breaks is bigger than the difference between limits.')
    }
  # Else, if breaks is a vector...
  } else {
    if (any(breaks < limits[1]) | any(breaks > limits[2])) {
      warning('Removing breaks outside the specified limits.')
      breaks <- breaks[which(breaks >= limits[1] & breaks <= limits[2])]
    }
  }
  # --- Code ---
  # If breaks is not a vector...
  if (length(breaks) == 1) {
    ### Set a new center = center - breaks/2. The new center has to be at a
    ### minimum distance of breaks/2 from the limits or be the upper limit
    ### itself.
    if (center < (limits[1] + (breaks/2))) {
      original.center <- center
      center <- limits[1] + (breaks/2)
    } else if (center > (limits[2] - (breaks/2))) {
      original.center <- center
      center <- limits[2]
    }
    center <- center - (breaks/2)
    ### Compute brk.low (from the lower limit to the new center, by breaks).
    if (limits[1] < center) {
      brk.low <- c(limits[1], seq(from = center, to = limits[1], by = -breaks))
      brk.low <- sort(unique(brk.low[which(brk.low >= limits[1])]),
                      decreasing = FALSE)
    } else {
      brk.low <- center
    }
    ### Compute brk.high (from the new center to the upper limit, by breaks).
    if (limits[2] > center) {
      brk.high <- c(limits[2], seq(from = center, to = limits[2], by = breaks))
      brk.high <- sort(unique(brk.high[which(brk.high <= limits[2])]),
                       decreasing = FALSE)
    } else {
      brk.high <- center
    }
    ### Pseudocenter: the new center + breaks/2.
    pseudo.center <- tail(brk.low, n = 1) + breaks/2
    ### Final breaks.
    final.breaks <- brk.labels <- sort(unique(c(brk.low, pseudo.center, brk.high)),
                                       decreasing = FALSE)
    ### Remove all labels but the limits and the pseudo.center.
    brk.labels[which(!(brk.labels %in% c(pseudo.center, limits)))] <- ""
    idx.pseudo.center <- which(brk.labels == pseudo.center)
    ### If the original.center was at a distance < breaks/2 from the limits,
    ### modify the labels and the pseudo.center.
    if (exists("original.center")) {
      brk.labels[idx.pseudo.center] <- pseudo.center <- original.center
    }
    ### If the pseudo.center == limit, remove the limits from the labels.
    if (pseudo.center %in% limits) brk.labels[-c(1, length(brk.labels))] <- ""
    ### If breaks is a vector...
  } else {
    ### Add limits to breaks.
    breaks <- sort(unique(c(limits, breaks)), decreasing = FALSE)
    ### The pseudocenter = original center.
    pseudo.center <- center
    ### Check which breaks element is the minimum value that is >= center.
    idx.bigger.than.center <- which(cumsum(breaks >= center) == 1)
    ### Set the new center to the previous element (or to the 1st element if
    ### idx.bigger.than.center == 1).
    if (idx.bigger.than.center > 1) center <- breaks[idx.bigger.than.center - 1]
    else center <- breaks[idx.bigger.than.center]
    ### brk.low (from the lower limit to the new center).
    brk.low <- breaks[1:which(breaks == center)]
    ### brk.high (from the new center + 1 to the upper limit).
    brk.high <- breaks[(which(breaks == center) + 1):length(breaks)]
    ### Final breaks and labels.
    final.breaks <- brk.labels <- sort(unique(c(brk.low, pseudo.center, brk.high)),
                                       decreasing = FALSE)
  }
  # Colours.
  # The colour of the center (last break of brk.low and first break of brk.high)
  # and the pseudo.center is the same (so these three values form a single
  # colour break).
  rampcol.mid <- rep(colorscale[3], times = 3)
  # If brk.low is more than just the center, get a different colour for each
  # break.
  if (length(brk.low) > 1) {
    rampcol.low <- colorRampPalette(colors = colorscale[1:2],
                                    space = "Lab")(length(brk.low)-1)
  } else rampcol.low <- character(0)
  # If brk.high is more than just the center and the limits[2], get a different
  # coloyr for each break.
  if (length(brk.high) > 2) {
    rampcol.high <- colorRampPalette(colors = colorscale[4:5],
                                     space = "Lab")(length(brk.high)-1)
  } else rampcol.high <- character(0)
  # Rampcolors is the vector with the final colours.
  rampcolors <- c(rampcol.low, rampcol.mid, rampcol.high)
  # Guide argument.
  guide <- ggplot2::guide_coloursteps(even.steps = FALSE, show.limits = FALSE,
                                      title = "Beyondcell")
  # Output: ggplot2 colour scale.
  out <- scale_colour_stepsn(colours = scales::alpha(rampcolors, alpha = alpha),
                             breaks = final.breaks, labels = brk.labels,
                             values = scales::rescale(final.breaks, to = c(0, 1)),
                             na.value = scales::alpha(na.value, alpha = alpha),
                             limits = limits, guide = guide)
  return(out)
}

#' @title Breaks a string into several lines
#' @description This function breaks a string \code{x} formed by elements
#' separated by \code{split} into lines of length \code{line.length}.
#' @name BreakString
#' @param x String to be breaked, formed by elements separated by \code{split}.
#' @param split Character that separates the elements of \code{x}.
#' @param line.length Length of the lines into which \code{x} will be breaked.
#' @return A string with the same contents that \code{x} breaked in lines of
#' length \code{line.length}.
#' @examples
#' @export

BreakString <- function(x, split = ", ", line.length = 50) {
  # --- Checks ---
  # Check x.
  if (length(x) != 1 | !is.character(x)) {
    stop('x must be a single string.')
  }
  # Check split.
  if (length(split) != 1 | !is.character(split)) {
    stop('split must be a single string.')
  }
  # Check line.length.
  if (length(line.length) != 1 | line.length[1] < 1 | line.length[1]%%1 != 0) {
    stop('line.length must be a single integer > 0.')
  }
  # --- Code ---
  if (nchar(x) <= line.length) final.x <- x
  else {
    split.x <- unlist(strsplit(x, split = split))
    # Length of each element in x + length(split).
    n.char <- sapply(split.x, function(y) nchar(y) + length(split))
    # Line in which each element of x will be printed.
    n.line <- cumsum(n.char)%/%line.length + 1
    # Separate each element within a line with split; and eachline with split + "\n".
    final.x <- paste0(sapply(1:max(n.line), function(i) {
      sub.x <- paste0(names(which(n.line == i)), collapse = split)
      return(sub.x)
    }), collapse = paste0(gsub(pattern = "^\\s+|\\s+$", replacement = "",
                               x = split), "\n")) # Trim
  }
  return(final.x)
}

#' @title Returns a dataframe with information about the input drugs
#' @description This function searches information about the inputted drugs in
#' the pre-loaded \code{\link[beyondcell]{beyondcell}} matrices and returns a
#' \code{data.frame} with their drug synonyms and MoAs.
#' @name FindDrugs
#' @param bc \code{\link[beyondcell]{beyondcell}} object.
#' @param x A character vector with drug names and/or sig IDs.
#' @details The output \code{data.frame} has the following columns:
#' \itemize{
#' \item{\code{Original_Name}}: Inputted drug name.
#' \item{\code{bc_Name}}: Drug name used in \code{bc}.
#' \item{\code{Preferred_Name}}: Drug name used in \code{beyondcell} plots.
#' \item{\code{Name}}: Other drug names.
#' \item{\code{sig_id}}: Signature ID.
#' \item{\code{Preferred_and_sig}}: \code{Preferred_Name} (or alternatively
#' \code{bc_Name}) and \code{sig_id}.
#' \item{\code{MoA}}: Mechanism(s) of action.
#' }
#' @return A \code{data.frame}.
#' @examples
#' @export

FindDrugs <- function(bc, x) {
  # --- Checks ---
  # Check that bc is a beyondcell object.
  if (class(bc) != "beyondcell") stop('bc must be a beyondcell object.')
  # Check x.
  if (!is.character(x)) stop('x must be a character vector.')
  # --- Code ---
  # bc signatures.
  sigs <- rownames(bc@normalized)
  # Match x with bc signatures and get the indices of matching elements.
  indices <- lapply(x, function(y) {
    idx <- match(toupper(y), table = toupper(sigs), nomatch = 0)
    if (idx == 0) {
      idx <- unique(match(drugInfo$sig_id[drugInfo$Name == toupper(y)],
                          table = sigs))
    }
    return(idx[!is.na(idx)])
  })
  # Original names (x) and bc names (sigs).
  df <- data.frame(Original_Name = unlist(sapply(seq_along(x), function(i) {
    rep(x[i], times = length(indices[[i]]))
  })), sig_id = unlist(sapply(indices, function(z) sigs[z])))
  # Get the names and pathways of the selected signatures.
  info <- subset(drugInfo, subset = sig_id %in% df$sig_id)
  info <- aggregate(.~ sig_id, data = info, na.action = NULL, FUN = function(y) {
    paste(na.omit(unique(y)), collapse = ", ")
  })
  # Merge df and info.
  df <- merge(df, info[, c("sig_id", "Name", "Preferred_Name", "MoA")],
              by = "sig_id", all.x = TRUE)
  # Add bc_Name column and remove names that are not sig IDs from sig_id column.
  df$bc_Name <- df$sig_id
  df$sig_id[!startsWith(df$sig_id, prefix = "sig_")] <- NA
  # Create Preferred_and_sig column: Preferred name and sig_id.
  df$Preferred_and_sig <- sapply(1:nrow(df), function(i) {
    name <- ifelse(test = !is.na(df$Preferred_Name[i]), yes = df$Preferred[i],
                   no = df$bc_Name[i])
    sig <- ifelse(test = !is.na(df$sig_id[i]),
                  yes = paste0(" (", df$sig_id[i], ")"), no = "")
    return(paste0(name, sig))
  })
  # Reorder df.
  df <- df[c("Original_Name", "bc_Name", "Preferred_Name", "Name",
             "sig_id", "Preferred_and_sig", "MoA")]
  return(df)
}
