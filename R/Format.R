#' @title Returns a vector of 5 colours
#' @description This function gets a colour palette and returns a vector with 5
#' colour values to be used in
#' \code{\link[beyondcell]{center_scale_colour_stepsn}}.
#' @name get_colour_steps
#' @importFrom viridis viridis
#' @importrom RColorBrewer brewer.pal
#' @param colorscale Either a \code{viridis}, \code{RColorBrewer} or a custom
#' palette of 3 colours (low, medium and high). If \code{colorscale = NULL}
#' (default), the function returns \code{beyondcell}'s own palette.
#' @return A vector with 5 colour values.
#' @examples
#' @export

get_colour_steps <- function(colorscale = NULL) {
  # --- Checks and Code ---
  default <- c("#1D61F2", "#83A8F7", "#F7F7F7", "#FF9CBB", "#DA0078")
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

#' @title Creates a centred sequential binned colour gradient
#' @description This function creates a sequential binned colour gradient
#' (low-mid-high) centred around \code{center}.
#' @name center_scale_colour_stepsn
#' @import ggplot2
#' @import scales
#' @param x A numeric vector. It can contain \code{NA}s.
#' @param colorscale A vector with 5 colours that can be obtained using
#' \code{\link[beyondcell]{get_colour_steps}}.
#' @param alpha Transparency level between 1 (not transparent) and 0 (fully
#' transparent).
#' @param na.value Colour to use for missing values.
#' @param limits Vector with the desired limits.
#' @param center A single number indicating the centre of the \code{colorscale}.
#' If \code{center = NULL} (default), the centre is set to the middle point of
#' \code{x}.
#' @param breaks A single number indicating the break size of the
#' \code{colorscale}. Alternatively, it can be a vector with the desired breaks
#' (which don't have to be symmetric or equally distributed).
#' @param aesthetics (\code{\link[ggplot2]{scale_colour_stepsn}}'s 
#' \code{aesthetics}) Character string or vector of character strings listing 
#' the name(s) of the aesthetic(s) that this scale works with. This can be 
#' useful, for example, to apply colour settings to the \code{colour} and 
#' \code{fill} aesthetics at the same time, via \code{aesthetics = c("colour", 
#' "fill")}.
#' @return A centred sequential binned colour gradient that can be used to
#' colour \code{\link[ggplot2]{ggplot2}} objects.
#' @examples
#' @export

center_scale_colour_stepsn <- function(x, colorscale, alpha = 0.7,
                                       na.value = "grey50", limits = c(NA, NA),
                                       center = NULL, breaks = 0.1, 
                                       aesthetics = "colour") {
  # --- Checks ---
  # Check x.
  if (!is.numeric(x)) {
    stop('x must be a numeric vector.')
  }
  range.values <- pretty(x)
  # Check colorscale.
  if (length(colorscale) != 5 |
      !tryCatch(is.matrix(col2rgb(colorscale)), error = function(cond) FALSE)) {
    stop('colorscale must contain exactly 5 colours.')
  }
  # Check alpha.
  if (length(alpha) != 1 | alpha[1] < 0 | alpha[1] > 1) {
    stop('alpha must be a positive number between 0 and 1.')
  }
  # Check na.value.
  if (!tryCatch(is.matrix(col2rgb(na.value)), error = function(cond) FALSE)) {
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
    ### Define center's limits at a maximum distance of breaks/2.
    limits.center <- c(max(limits[1], center - (breaks/2)), 
                       min(limits[2], center + (breaks/2)))
    ### Compute brk.low (from the global lower limit to the center's lower 
    ### limit, by breaks).
    if (limits.center[1] == limits[1]) brk.low <-limits[1]
    else {
      brk.low <- c(limits[1], 
                   seq(from = limits.center[1], to = limits[1], by = -breaks))
      brk.low <- sort(unique(brk.low[which(brk.low >= limits[1])]),
                      decreasing = FALSE)
    }
    ### Compute brk.high (from the center's upper limit to the global upper 
    ### limit, by breaks).
    if (limits.center[2] == limits[2]) brk.high <- limits[2]
    else {
      brk.high <- c(limits[2], 
                    seq(from = limits.center[2], to = limits[2], by = breaks))
      brk.high <- sort(unique(brk.high[which(brk.high <= limits[2])]),
                       decreasing = FALSE)
    }
    ### Final breaks.
    final.breaks <- 
      brk.labels <- sort(unique(na.omit(c(brk.low, center, brk.high))),
                         decreasing = FALSE)
    ### Remove all labels but the limits and the center.
    idx.center <- which(brk.labels == center)
    idx.limits <- which(brk.labels %in% limits)
    brk.labels[-c(idx.center, idx.limits)] <- ""
    ### If the center is too close to one limit but it is not the limit, don't 
    ### print its label.
    if ((limits.center[1] == limits[1] | limits.center[2] == limits[2]) &
        !center %in% limits.center) brk.labels[idx.center] <- ""
    ### If breaks is a vector...
  } else {
    ### Add limits to breaks.
    breaks <- sort(unique(c(limits, breaks)), decreasing = FALSE)
    ### Define center's limits (closer upper and lower breaks).
    idx.lower.than.center <- max(which(breaks <= center))
    idx.bigger.than.center <- min(which(breaks >= center))
    ### Compute brk.low (from the global lower limit to the center's lower 
    ### limit.
    brk.low <- breaks[1:idx.lower.than.center]
    ### Compute brk.high (from the center's upper limit to the global upper 
    ### limit.
    brk.high <- breaks[idx.bigger.than.center:length(breaks)]
    ### Final breaks and labels.
    final.breaks <- brk.labels <- sort(unique(c(brk.low, center, brk.high)),
                                       decreasing = FALSE)
  }
  # Colours.
  # The colour of the center and its limits is the same (so these three values 
  # form a single colour break).
  rampcol.mid <- rep(colorscale[3], times = 3)
  # If brk.low is more than just the lower limit, get a different colour for 
  # each break.
  if (length(brk.low) > 1) {
    rampcol.low <- colorRampPalette(colors = colorscale[1:2],
                                    space = "Lab")(length(brk.low)-1)
  } else rampcol.low <- character(0)
  # If brk.high is  more than just the upper limit, get a different colour for 
  # each break.
  if (length(brk.high) > 1) {
    rampcol.high <- colorRampPalette(colors = colorscale[4:5],
                                     space = "Lab")(length(brk.high)-1)
  } else rampcol.high <- character(0)
  # Rampcolors is the vector with the final colours.
  rampcolors <- c(rampcol.low, rampcol.mid, rampcol.high)
  # Guide argument.
  guide <- ggplot2::guide_coloursteps(even.steps = FALSE, show.limits = FALSE,
                                      title = "BCS", label.vjust = 0.5)
  # Output: ggplot2 colour scale.
  out <- scale_colour_stepsn(colours = scales::alpha(rampcolors, alpha = alpha),
                             breaks = final.breaks, labels = brk.labels,
                             values = scales::rescale(final.breaks, to = c(0, 1)),
                             na.value = scales::alpha(na.value, alpha = alpha),
                             limits = limits, guide = guide, 
                             aesthetics = aesthetics)
  return(out)
}

#' @title Breaks a string into several lines
#' @description This function breaks a string \code{x} formed by elements
#' separated by \code{split} into lines of length \code{line.length}.
#' @name BreakString
#' @param x String to be broken, formed by elements separated by \code{split}.
#' @param split Character that separates the elements of \code{x}.
#' @param line.length Length of the lines into which \code{x} will be broken.
#' @return A string with the same content as \code{x} broken in lines of
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
