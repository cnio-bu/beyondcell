#' @title Subsets a beyondcell object
#' @description This function subsets a \code{\link[beyondcell]{beyondcell}}
#' object based on signature names, cells and/or the maximum proportion of
#' \code{NaN} values desired in each signature and/or cell. See Details for more
#' information.
#' @name bcSubset
#' @param bc \code{beyondcell} object.
#' @param signatures Vector with the names of the signatures to subset by. If
#' \code{signatures = NULL}, all signatures will be kept.
#' @param bg.signatures Vector with the names of the background signatures to
#' subset by. If \code{bg.signatures = NULL}, all background signatures will be
#' kept.
#' @param cells Vector with the names of the cells to subset by. If
#' \code{cells = NULL}, all cells will be kept.
#' @param nan.sigs Maximum proportion of \code{NaN} values per signature in the
#' output \code{beyondcell} object. All signatures with a proportion of
#' \code{NaN > nan.sig} will be removed.
#' @param nan.cells Maximum proportion of \code{NaN} values per cell in the
#' output \code{beyondcell} object. All cells with a proportion of
#' \code{NaN > nan.cells} will be removed.
#' @details This function can subset a \code{beyondcell} object using its 5
#' parameters alone or in combination.
#'
#' So, for example, if you specify \code{signatures} and \code{cells}, the
#' resulting \code{beyondcell} object (except the \code{@@background} slot) will
#' be subsetted according to those vectors. The slot \code{@@background} will be
#' only subsetted according to \code{cells}. If you want to subset it by
#' signatures as well, you must specify a value for \code{bg.signatures}.
#'
#' On the other hand, if you specify \code{cells} and \code{nan.sigs}, the
#' output \code{beyondcell} object will keep the selected cells and those
#' signatures with a proportion of \code{NaN} values below or equal to
#' \code{nan.sigs}. Note that \code{nan.sigs} and \code{nan.cells} arguments
#' also subset the signatures and cells that meet the criteria in
#' \code{@@background} slot.
#'
#' Finally, if you specify all parameters, the result will keep those signatures
#' and cells of interest with a proportion of \code{NaN} below or equal to
#' \code{nan.sigs} and \code{nan.cells}, respectively.
#' @return A subsetted \code{beyondcell} object.
#' @examples
#' @export

bcSubset <- function(bc, signatures = NULL, bg.signatures = NULL, cells = NULL,
                     nan.sigs = 1, nan.cells = 1) {
  # --- Checks ---
  # Check that bc is a beyondcell object.
  if (class(bc) != "beyondcell") stop('bc must be a beyondcell object.')
  # Check signatures.
  if (is.null(signatures)) {
    pass.sigs <- rownames(bc@scaled)
  } else {
    in.signatures <- signatures %in% rownames(bc@scaled)
    if (all(!in.signatures)) {
      stop('None of the specified signatures were found.')
    } else if (any(!in.signatures)) {
      warning(paste0('These signatures were not found in bc: ',
                     paste0(unique(signatures[!in.signatures]),
                            collapse = ", "), '.'))
    }
    pass.sigs <- unique(signatures[in.signatures])
  }
  # Check bg.signatures
  if (is.null(bg.signatures)) {
    pass.sigs.bg <- rownames(bc@background)
  } else {
    in.signatures.bg <- bg.signatures %in% rownames(bc@background)
    if (all(!in.signatures.bg)) {
      stop('None of the specified bg.signatures were found.')
    } else if (any(!in.signatures.bg)) {
      warning(paste0('These bg.signatures were not found in bc: ',
                     paste0(unique(bg.signatures[!in.signatures.bg]),
                            collapse = ", "), '.'))
    }
    pass.sigs.bg <- unique(bg.signatures[in.signatures.bg])
  }
  # Check cells.
  if (is.null(cells)) {
    pass.cells <- colnames(bc@scaled)
  } else {
    in.cells <- cells %in% colnames(bc@scaled)
    if (all(!in.cells)) {
      stop('None of the specified cells were found.')
    } else if (any(!in.cells)) {
      warning(paste0('These cells were not found in bc matrices: ',
                     paste0(unique(cells[!in.cells]), collapse = ", "), '.'))
    }
    pass.cells <- unique(cells[in.cells])
  }
  # If there is spatial info, check cells in the slices too.
  is.spatial <- exists(x = "images", where = bc@SeuratInfo)
  if (is.spatial) {
    if(length(bc@SeuratInfo$images) > 0) {
    spots <- lapply(bc@SeuratInfo$images, FUN = function(x) {
      return(rownames(x@coordinates))
    })
    spots <- unique(unlist(spots))
    in.pass.cells <- pass.cells %in% spots
    if (all(!in.pass.cells)) {
      stop('None of the specified cells were found in the image slices.')
    } else if (any(!in.pass.cells)) {
      warning(paste0('These cells were not found in the image slices: ',
                     paste0(pass.cells[!in.pass.cells], collapse = ", "), '.'))
      }
    }
  }
  # Check nan.sigs.
  if (length(nan.sigs) != 1 | nan.sigs[1] < 0 | nan.sigs[1] > 1) {
    stop('nan.sigs must be a positive number between 0 and 1.')
  }
  pass.nan.sigs <- rownames(bc@scaled)[apply(bc@scaled, 1, function(x) {
    sum(is.na(x)) <= ncol(bc@scaled) * nan.sigs
  })]
  pass.nan.sigs.bg <- rownames(bc@background)[apply(bc@background, 1,
    function(y) sum(is.na(y)) <= ncol(bc@background) * nan.sigs)]
  # Check nan.cells.
  if (length(nan.cells) != 1 | nan.cells[1] < 0 | nan.cells[1] > 1) {
    stop('nan.cells must be a positive number between 0 and 1.')
  }
  pass.nan.cells <- colnames(bc@scaled)[apply(bc@scaled, 2, function(y) {
    sum(is.na(y)) <= nrow(bc@scaled) * nan.cells
  })]
  # Check regression and subset order.
  reg.order <- bc@regression$order
  reg.order.bg <- bc@regression$order.background
  if (any(!(reg.order %in% c("regression", "subset", ""))) |
      reg.order[1] == "" & reg.order[2] != "" | length(reg.order) != 2 |
      (identical(reg.order[1], reg.order[2]) & reg.order[1] != "") |
      (!identical(reg.order.bg, c("", "")) & !identical(reg.order, reg.order.bg))) {
    warning(paste('Corrupt beyondcell object. Restoring original object before',
                  'subsetting...'))
    bc <- bcRecompute(bc, slot = "data")
    bc@background <- matrix(ncol = 0, nrow = 0)
    bc@regression <- list(order = c("subset", ""), vars = NULL,
                          order.background = rep("", times = 2))
    reg.order <- rep("", 2)
  } else {
    ### If bc was previously subsetted and then regressed, raise an error.
    if (identical(reg.order, c("subset", "regression"))) {
      stop(paste('bc was previously subsetted and then regressed. Please run',
                 'bcSubset() on a new beyondcell object created with',
                 'CreatebcObject(bc).'))
      ### Else if reg.order == c("", "") or reg.order == c("regression", ""), add
      ### "subset" to bc@regression$order.
    } else if (identical(reg.order, rep("", 2)) |
               identical(reg.order, c("regression", ""))) {
      bc@regression$order[match("", table = reg.order)] <- "subset"
      ### Else if the last step in reg.order is "subset", raise a warning.
    } else if (tail(reg.order[which(reg.order != "")], n = 1) == "subset") {
      warning('bc is an already subsetted object.')
    }
  }
  # --- Code ---
  # Intersect pass.sig and pass.nan.sig (signatures that pass both filters).
  final.sigs <- intersect(pass.sigs, pass.nan.sigs)
  # Intersect pass.sig.bg and pass.nan.sig.bg (bg.signatures that pass both
  # filters).
  final.sigs.bg <- intersect(pass.sigs.bg, pass.nan.sigs.bg)
  # Intersect pass.cells and pass.nan.cells (cells that pass both filters).
  final.cells <- intersect(pass.cells, pass.nan.cells)
  # If final.sigs == rownames(bc@scaled),
  # final.sigs.bg == rownames(bc@background) and
  # final.cells == colnames(bc@scaled) == colnames(bc@background),
  # the output is equal to the input. Else, subset.
  if (!identical(sort(final.sigs, decreasing = FALSE),
                 sort(rownames(bc@normalized), decreasing = FALSE)) |
      !identical(sort(final.sigs.bg, decreasing = FALSE),
                 sort(rownames(bc@background), decreasing = FALSE)) |
      !identical(sort(final.cells, decreasing = FALSE),
                 sort(colnames(bc@normalized), decreasing = FALSE)) |
      !identical(sort(final.cells, decreasing = FALSE),
                 sort(colnames(bc@background), decreasing = FALSE))) {
    if (length(final.sigs) == 0) stop("No signature met all filter criteria.")
    else if (length(final.cells) == 0) stop("No cell met all filter criteria.")
    else {
      bc@normalized <- round(bc@normalized[final.sigs, final.cells,
                                           drop = FALSE], digits = 2)
      bc <- bcRecompute(bc, slot = "normalized")

      ### If there is spatial info, then subset the cells in the slices too
      if (is.spatial) {
        bc@SeuratInfo$images <- lapply(bc@SeuratInfo$images, FUN = function(x) {
          slice.cells <- rownames(x@coordinates)
          final.slice.cells <- intersect(slice.cells, final.cells)
          x@coordinates <- x@coordinates[final.slice.cells, , drop = FALSE]
          return(x)
        })
      }
    }
    if (any(dim(bc@background) != 0)) {
      if (length(final.sigs.bg) == 0) {
        stop("No background signature met all filter criteria.")
      } else {
        bc@background <- bc@background[final.sigs.bg, final.cells, drop = FALSE]
        bc@regression$order.background <- bc@regression$order
      }
    }
  } else {
    warning('bc was not subsetted.')
    bc@regression$order <- reg.order
  }
  return(bc)
}

#' @title Regresses out unwanted effects from BCS
#' @description This function regresses out unwanted effects from normalized
#' beyondcell scores (BCS).
#' @name bcRegressOut
#' @importFrom DMwR knnImputation
#' @param bc \code{\link[beyondcell]{beyondcell}} object.
#' @param vars.to.regress Vector of metadata columns to regress out the BCS.
#' @param k.neighbors (\code{\link[DMwR]{knnImputation}}'s \code{k}) Number of 
#' nearest neighbors to use.
#' @param add.DSS Use background BCS computed with \code{DSS} signatures
#' (\code{add.DSS = TRUE}) or just use the signatures included in the \code{bc}
#' object (\code{add.DSS = FALSE}) to do the imputation of \code{NaN} BCS. If 
#' the number of drugs in \code{bc} (excluding pathways) is <= 20, it is 
#' recomended to set \code{add.DSS = TRUE}. Note that if \code{add.DSS = TRUE}, 
#' the regression and subset steps that have been applied on \code{bc} will 
#' also be applied on the background BCS.
#' @return Returns a \code{beyondcell} object with regressed normalized BCS,
#' regressed scaled BCS and regressed switch points.
#' @examples
#' @export

bcRegressOut <- function(bc, vars.to.regress, k.neighbors = 10,
                         add.DSS = FALSE) {
  # --- Checks ---
  # Check that bc is a beyondcell object.
  if (class(bc) != "beyondcell") stop('bc must be a beyondcell object.')
  # Check vars.to.regress.
  in.vars.to.regress <- !is.null(vars.to.regress) &
    vars.to.regress %in% colnames(bc@meta.data)
  if (all(!in.vars.to.regress)) {
    stop('vars.to.regress not found.')
  } else if (any(!in.vars.to.regress)) {
    stop(paste0('Some vars.to.regress not found: ',
                paste0(vars.to.regress[!in.vars.to.regress], collapse = ", "),
                '.'))
  }
  vars <- unique(vars.to.regress[in.vars.to.regress])
  # Check regression and subset order.
  reg.order <- bc@regression$order
  reg.order.bg <- bc@regression$order.background
  reg.vars <- bc@regression$vars
  if (any(!(reg.order %in% c("regression", "subset", ""))) |
      reg.order[1] == "" & reg.order[2] != "" | length(reg.order) != 2 |
      (identical(reg.order[1], reg.order[2]) & reg.order[1] != "") |
      (!identical(reg.order.bg, c("", "")) & !identical(reg.order, reg.order.bg))) {
    warning(paste('Corrupt beyondcell object. Restoring original object before',
                  'regressing...'))
    bc <- suppressMessages(bcRecompute(bc, slot = "data"))
    bc@background <- matrix(ncol = 0, nrow = 0)
    reg.order <- reg.order.bg <- rep("", 2)
    reg.vars <- NULL
  } else {
    ### If bc was previously regressed and then subsetted, raise an error.
    if (identical(reg.order, c("regression", "subset"))) {
      stop(paste('bc was previously regressed and then subsetted. Please run',
                 'bcSubset() on a new beyondcell object created with',
                 'CreatebcObject(bc).'))
    ### Else if the last step in reg.order is "regression", raise a warning.
    } else if (!identical(reg.order, rep("", 2))) {
      if (tail(reg.order[which(reg.order != "")], n = 1) == "regression") {
        warning('bc is an already regressed object.')
        vars <- unique(c(vars, reg.vars))
        bc <- suppressMessages(bcRecompute(bc, slot = "data"))
        bc@regression <- list(order = rep("", 2),  vars = NULL,
                              order.background = rep("", 2))
        if (any(dim(bc@background) != 0)) {
          message('Restoring pre-regressed background matrix...')
          gs.background <- suppressMessages(
            GetCollection(DSS, n.genes = bc@n.genes, mode = bc@mode,
                          include.pathways = FALSE))
          background <- suppressWarnings(suppressMessages(
            bcScore(bc@expr.matrix, gs = gs.background, expr.thres = bc@thres)))
          bc@background <- background@normalized
        }
        if ("subset" %in% reg.order) {
          bc <- suppressWarnings(
            bcSubset(bc, signatures = rownames(bc@normalized), 
                     bg.signatures = rownames(bc@background),
                     cells = colnames(bc@normalized)))
        }
        reg.order[reg.order == "regression"] <- ""
        reg.order.bg[reg.order.bg == "regression"] <- ""
      }
    }
  }
  bc@regression <- list(order = reg.order, vars = reg.vars, 
                        order.background = reg.order.bg)
  # Check k.neighbors.
  if (!is.numeric(k.neighbors)) stop('k.neighbors must be numeric.')
  if (length(k.neighbors) != 1 | k.neighbors[1]%%1 != 0 | k.neighbors[1] < 1) {
    stop('k.neighbors must be a positive integer.')
  }
  # Check add.DSS.
  cells <- colnames(bc@normalized)
  sigs <- rownames(bc@normalized)
  not.paths <- !(sigs %in% names(pathways))
  drugs <- sigs[not.paths]
  n.drugs <- sum(not.paths)
  n.complete.normalized <- sum(complete.cases(t(bc@normalized[drugs, , 
                                                              drop = FALSE])))
  is.complete.normalized <- n.complete.normalized == length(cells)
  if (length(add.DSS) != 1 | !is.logical(add.DSS)) {
    stop('add.DSS must be TRUE or FALSE.')
  } else if(add.DSS & is.complete.normalized) {
    warning('No NaN values were found in bc@normalized. add.DSS is deprecated.')
    add.DSS <- FALSE
  } else if (!add.DSS) {
    if (!is.complete.normalized) {
      if (n.drugs <= 10) {
        stop(paste('Only', n.drugs, 'drug signatures (excluding pathways) are',
                   'present in the bc object, please set add.DSS = TRUE.'))
      } else if (n.drugs <= 20) {
        warning(paste('Computing an UMAP reduction for', n.drugs,
                      'drugs. We recommend to set add.DSS = TRUE when the', 
                      'number of signatures (excluding pathways) is below or', 
                      'equal to 20.'))
      }
    }
    ### Complete cases for normalized BCS.
    if (k.neighbors >= n.complete.normalized) {
      stop(paste0('k.neighbors must be lower than the number of complete ', 
                  'cases in @normalized slot: ', n.complete.normalized, 
                  '.'))
    }
    ### Complete cases for background BCS.
    if (all(dim(bc@background) == 0)) n.complete.bg <- length(cells)
    else n.complete.bg <- sum(complete.cases(t(bc@background)))
    if (k.neighbors >= n.complete.bg) {
      stop(paste0('k.neighbors must be lower than the number of complete ', 
                  'cases in @background slot: ', n.complete.bg, '.'))
    }
  }
  # --- Code ---
  if (add.DSS) {
    ### DSS (background) BCS.
    if (!identical(sort(rownames(bc@background), decreasing = FALSE),
                   sort(unique(DSS@info$IDs), decreasing = FALSE)) |
        !identical(sort(colnames(bc@background), decreasing = FALSE),
                   sort(cells, decreasing = FALSE)) |
        !identical(bc@regression$order, bc@regression$order.background)) {
      message('Computing background BCS using DSS signatures...')
      ### Genesets.
      gs.background <- suppressMessages(
        GetCollection(DSS, n.genes = bc@n.genes, mode = bc@mode, 
                      include.pathways = FALSE))
      ### BCS.
      background <- suppressWarnings(suppressMessages(
        bcScore(bc@expr.matrix, gs = gs.background, expr.thres = bc@thres)))
      ### Add metadata.
      background@meta.data <- background@meta.data[, -c(1:ncol(background@meta.data)), 
                                                   drop = FALSE]
      background <- bcAddMetadata(background, metadata = bc@meta.data)
      ### Subset if needed.
      if (reg.order[1] == "subset") {
        background <- bcSubset(background, cells = cells)
      }
      ### Add background@normalized to bc@background.
      bc@background <- background@normalized
      ### Update reg.order.bg.
      reg.order.bg <- reg.order
    } else {
      message('Background BCS already computed. Skipping this step.')
    }
    ### Add background to bc.
    all.rows <- unique(c(drugs, rownames(bc@background)))
    merged.score <- rbind(bc@normalized, 
                          bc@background[, cells, drop = FALSE])[all.rows, , 
                                                                drop = FALSE]
    bc.merged <- beyondcell(normalized = merged.score)
    ### Complete cases for merged BCS.
    n.complete.merged <- sum(complete.cases(t(bc.merged@normalized)))
    if (k.neighbors >= n.complete.merged) {
      stop(paste0('k.neighbors must be lower than the total number of ', 
                  'complete cases in @normalized and @background slots: ', 
                  n.complete.merged, '.'))
    }
  } else {
    ### No background BCS.
    message(paste('DSS background not computed. The imputation will be', 
                  'computed with just the drugs (not pathways) in the', 
                  'beyondcell object.'))
    bc.merged <- beyondcell(normalized = bc@normalized[drugs, , drop = FALSE])
  }
  # Latent data.
  latent.data <- bc@meta.data[cells, vars, drop = FALSE]
  # Impute normalized BCS matrix if necessary
  if (!is.complete.normalized) {
    message('Imputing normalized BCS...')
    result <- t(DMwR::knnImputation(t(bc.merged@normalized), k = k.neighbors, 
                                    scale = FALSE, meth = "weighAvg"))
  } else {
    message('No imputation needed for bc@normalized.')
    result <- bc.merged@normalized
  }
  imputation <- result[drugs, cells, drop = FALSE]
  # Limma formula.
  fmla <- as.formula(object = paste('bcscore ~', paste(vars, collapse = '+')))
  # Compute regression and save it in bc@normalized.
  message('Regressing scores...')
  total <- nrow(imputation)
  pb <- txtProgressBar(min = 0, max = total, style = 3, file = stderr())
  bins <- ceiling(total / 100)
  normalized.regressed <- t(apply(cbind(seq_len(nrow(imputation)),
                                        imputation), 1, function(x) {
                            regression.mat <- cbind(latent.data, x[-1])
                            colnames(regression.mat) <- c(colnames(latent.data), "bcscore")
                            qr <- lm(fmla, data = regression.mat, qr = TRUE, na.action = na.exclude)$qr
                            resid <- qr.resid(qr = qr, y = x[-1])
                            ### Update the progress bar.
                            if (x[1]%%bins == 0 | x[1] == total) {
                              Sys.sleep(0.1)
                              setTxtProgressBar(pb, value = x[1])
                            }
                          ### Return residues.
                          return(resid)
                        }))
  bc@normalized <- round(rbind(bc@normalized[!not.paths, , drop = FALSE], 
                               normalized.regressed)[sigs, cells, drop = FALSE], 
                         digits = 2)
  # Close the progress bar.
  Sys.sleep(0.1)
  close(pb)
  # Recompute the beyondcell object
  bc <- bcRecompute(bc, slot = "normalized")
  # Add "regression" step to bc@regression$order.
  reg.order[grep("^$", reg.order)[1]] <- "regression"
  bc@regression$order <- reg.order
  # Add vars.to.regress to bc@regression$vars.
  bc@regression$vars <- vars
  # Regress the background, if needed.
  if (any(dim(bc@background) != 0)) {
    if (!add.DSS) {
      is.complete.bg <- all(complete.cases(t(bc@background)))
      if (!is.complete.bg) {
        message('Imputing background BCS...')
        imputation.bg <- t(DMwR::knnImputation(t(bc@background), k = k.neighbors,
                                               scale = FALSE, meth = "weighAvg"))
      } else {
        message('No imputation needed for bc@background.')
        imputation.bg <- bc@background
      }
    } else {
      message('Background BCS already imputed.')
      imputation.bg <- result[unique(DSS@info$IDs), cells, drop = FALSE]
    }
    message('Regressing background BCS...')
    total.bg <- nrow(imputation.bg)
    pb.bg <- txtProgressBar(min = 0, max = total.bg, style = 3, file = stderr())
    bins.bg <- ceiling(total.bg / 100)
    background.regressed <- t(apply(cbind(seq_len(nrow(imputation.bg)),
                                          imputation.bg), 1, function(y) {
                                regression.mat <- cbind(latent.data, y[-1])
                                colnames(regression.mat) <- c(colnames(latent.data), "bcscore")
                                qr <- lm(fmla, data = regression.mat, qr = TRUE, na.action = na.exclude)$qr
                                resid <- qr.resid(qr = qr, y = y[-1])
                                ### Update the progress bar.
                                if (y[1]%%bins == 0 | y[1] == total.bg) {
                                  Sys.sleep(0.1)
                                  setTxtProgressBar(pb.bg, value = y[1])
                                }
                              ### Return residues.
                              return(resid)
                              }))
    bc@background <- round(background.regressed[, cells, drop = FALSE], 
                           digits = 2)
    bc@regression$order.background <- bc@regression$order
    # Close the background progress bar.
    Sys.sleep(0.1)
    close(pb.bg)
    # Add "regression" step to bc@regression$order.background.
    reg.order.bg[grep("^$", reg.order.bg)[1]] <- "regression" 
    bc@regression$order.background <- reg.order.bg
  }
  # Output.
  return(bc)
}

#' @title Recomputes a beyondcell object
#' @description This function recomputes a \code{\link[beyondcell]{beyondcell}}
#' object using the matrix stored in the slot \code{@@data} (original scores) or
#' \code{@@normalized} (which can contain subsetted and/or regressed scores).
#' Columns added with \code{\link[beyondcell]{bcAddMetadata}} are preserved,
#' except if they define therapeutic clusters. Important: \code{bc@background}
#' remains the same, while \code{bc@ranks} and \code{bc@reductions} are removed.
#' @name bcRecompute
#' @import scales
#' @param bc \code{beyondcell} object.
#' @param slot Score matrix to recompute the \code{beyondcell} object. Either
#' \code{"data"} or \code{"normalized"}.
#' @return A recomputed \code{beyondcell} object.
#' @examples
#' @export

bcRecompute <- function(bc, slot = "data") {
  # --- Checks ---
  # Check that bc is a beyondcell object.
  if (class(bc) != "beyondcell") stop('bc must be a beyondcell object.')
  # Check slot.
  if (length(slot) != 1 | !(slot[1] %in% c("data", "normalized"))) {
    stop('slot must be "data" or "normalized".')
  }
  # --- Code ---
  if (slot == "data") {
    # bc@normalized = bc@data.
    bc@normalized <- bc@data
  } else if (slot == "normalized") {
    bc@data <- bc@normalized <- round(bc@normalized, digits = 2)
  }
  # Recompute scaled BCS.
  scaled <- round(t(apply(bc@normalized, 1, scales::rescale, to = c(0, 1))),
                  digits = 2)
  rownames(scaled) <- rownames(bc@normalized)
  colnames(scaled) <- colnames(bc@normalized)
  bc@scaled <- scaled
  # Recompute switch points.
  bc@switch.point <- SwitchPoint(bc)
  # Remove ranks and reductions.
  if (!is.null(names(bc@ranks))) message('Removing @ranks slot...')
  if (!is.null(names(bc@reductions))) message('Removing @reductions slot...')
  bc@ranks <- bc@reductions <- vector(mode = "list", length = 0)
  # Remove therapeutic clusters from bc@meta.data.
  therapeutic.clusters <- grep(pattern = "bc_clusters_res.",
                               x = colnames(bc@meta.data))
  if (length(therapeutic.clusters) > 0) {
    message('Removing therapeutic clusters...')
    bc@meta.data <- bc@meta.data[, -c(therapeutic.clusters), drop = FALSE]
  }
  bc@meta.data <- bc@meta.data[colnames(bc@normalized), , drop = FALSE]
  # Return bc object.
  return(bc)
}

#' @title Add new metadata to an existing beyondcell object
#' @description This function adds new metadata to an existing
#' \code{\link[beyondcell]{beyondcell}} object.
#' @name bcAddMetadata
#' @param bc \code{beyondcell} object.
#' @param metadata Matrix or dataframe with metadata to add. Rownames should be
#' cell names and colnames should be numeric or categorical variables.
#' @return Returns a \code{beyondcell} object with updated metadata.
#' @examples
#' @export

bcAddMetadata <- function(bc, metadata) {
  # --- Checks ---
  # Check that bc is a beyondcell object.
  if (class(bc) != "beyondcell") stop('bc must be a beyondcell object.')
  # Check metadata.
  if(!is.matrix(metadata) & !is.data.frame(metadata)) {
    stop('metadata must be must be a matrix or a data.frame.')
  }
  # Check that the rownames of bc@meta.data and the new metadata are the same.
  if (!identical(sort(rownames(bc@meta.data), decreasing = FALSE),
                 sort(rownames(metadata), decreasing = FALSE))) {
    stop('metadata and bc@meta.data rownames are not the same.')
  }
  # Check that columns in metadata are different from the existing columns in
  # bc@meta.data.
  in.metadata <- colnames(metadata) %in% colnames(bc@meta.data)
  if (any(in.metadata)) {
    warning(paste0('Some metadata columns are already present in bc@meta.data ', 
                   'slot and will be overwritten: ',
                   paste0(colnames(metadata)[in.metadata], collapse = ", "), 
                   '.'))
  }
  # --- Code ---
  metadata <- metadata[rownames(bc@meta.data), , drop = FALSE]
  bc@meta.data <- cbind(bc@meta.data, metadata)
  return(bc)
}

#' @title Merges two beyondcell objects
#' @description This function merges two \code{\link[beyondcell]{beyondcell}}
#' objects obtained from the same single-cell matrix using the same
#' \code{expr.thres} (see \code{\link[beyondcell]{bcScore}} for more
#' information). It binds signatures, not cells. If \code{bc1} is a subset of 
#' \code{bc2}, the resulting \code{beyondcell} object will contain just the 
#' subsetted cells and the \code{SeuratInfo} stored in \code{bc1}.
#' @name bcMerge
#' @importFrom plyr join
#' @param bc1 First \code{beyondcell} object to merge. Can be a subset of 
#' \code{bc2}.
#' @param bc2 Second \code{beyondcell} object to merge.
#' @param keep.bc.clusters Whether to keep \code{bc1} reductions or not. 
#' Choose \code{keep.bc.clusters = TRUE} when \code{bc2} is formed by functional
#' pathways and \code{bc1} is a drug-containing \code{beyondcell} object.
#' @return A merged \code{beyondcell} object.
#' @examples
#' @export

bcMerge <- function(bc1, bc2, keep.bc.clusters = TRUE) {
  # --- Checks ---
  # Check that bc1 and bc2 are beyondcell objects.
  if (class(bc1) != "beyondcell") stop('bc1 must be a beyondcell object.')
  if (class(bc2) != "beyondcell") stop('bc2 must be a beyondcell object.')
  # Check keep.bc.clusters.
  if (length(keep.bc.clusters) != 1 | !is.logical(keep.bc.clusters[1])) {
    stop('keep.bc.clusters must be TRUE or FALSE.')
  }
  # Check both thres.
  if (!identical(bc1@thres, bc2@thres)) {
    stop('bc objects weren\'t obtained using the same expression threshold.')
  }
  # Check both Seurat experiments.
  if (!identical(bc1@expr.matrix, bc2@expr.matrix)) {
    stop('bc objects weren\'t obtained from the same single-cell experiment.')
  }
  # Check subsetted cells.
  common.cells <- intersect(colnames(bc1@data), colnames(bc2@data))
  if (length(common.cells) == 0) {
    stop('bc1 and bc2 do not contain the same cells.')
  } else if (!all(common.cells %in% colnames(bc1@data))) {
    stop('bc1 is not a subset of bc2.')
  }
  # Check for duplicated signatures.
  duplicated.sigs <- intersect(rownames(bc1@data), rownames(bc2@data))
  if (length(duplicated.sigs) > 0) {
    identical.sigs <- sapply(duplicated.sigs, function(x) {
      identical(bc1@data[x, common.cells, drop = FALSE], 
                bc2@data[x, common.cells, drop = FALSE])
    })
    if (any(!identical.sigs)) {
      stop(paste0('Duplicated signatures: ',
                  paste0(duplicated.sigs[!identical.sigs], collapse = ", "),
                  ' without matching BCS in slot @data.'))
    }
  }
  # Check regression steps.
  if (!identical(bc1@regression, bc2@regression) & 
      (!identical(bc1@regression$order, c("subset", "")) & 
       !identical(bc2@regression$order, rep("", 2)))) {
    stop(paste('The two objects must be subsetted and/or regressed in the same',
               'order and with the same variales. Alternatively, bc1 can be a', 
               'subset of bc2 (with no regression).'))
  }
  # --- Code ---
  # Suset bc2
  bc2 <- suppressWarnings(suppressMessages(bcSubset(bc2, cells = common.cells)))
  # Create a new beyondcell object.
  bc <- beyondcell(expr.matrix = bc1@expr.matrix, SeuratInfo = bc1@SeuratInfo,
                   regression = bc1@regression,
                   n.genes = unique(c(bc1@n.genes, bc2@n.genes)),
                   mode = unique(c(bc1@mode, bc2@mode)), thres = bc1@thres)
  # rbind scaled BCS.
  bc@scaled <- unique(rbind(bc1@scaled[, common.cells, drop = FALSE], 
                            bc2@scaled[, common.cells, drop = FALSE]))
  # rbind normalized BCS.
  bc@normalized <- unique(rbind(bc1@normalized[, common.cells, drop = FALSE], 
                                bc2@normalized[, common.cells, drop = FALSE]))
  # rbind data.
  bc@data <- unique(rbind(bc1@data[, common.cells, drop = FALSE], 
                          bc2@data[, common.cells, drop = FALSE]))
  # Merge switch.points.
  bc@switch.point <- c(bc1@switch.point, bc2@switch.point)[rownames(bc@scaled)]
  # Merge meta.data.
  bc@meta.data <- suppressMessages(plyr::join(bc1@meta.data[common.cells, , drop = FALSE], 
                                              bc2@meta.data[common.cells, , drop = FALSE]))
  rownames(bc@meta.data) <- common.cells
  # If keep.bc.clusters, keep bc1 reductions.
  if (keep.bc.clusters) {
    bc@reductions <- bc1@reductions
  # Else, remove the therapeutic clusters from bc@meta.data.
  } else {
    therapeutic.clusters <- grep(pattern = "bc_clusters_res.", x = colnames(bc@meta.data))
    if (length(therapeutic.clusters) > 0) {
      bc@meta.data <- bc@meta.data[, -c(therapeutic.clusters), drop = FALSE]
    }
  }
  # Merge backgrounds.
  bg <- list(bc1 = as.data.frame(bc1@background), bc2 = as.data.frame(bc2@background))
  is.empty.bg <- sapply(bg, FUN = function(x) dim(x)[2] == 0)
  if (all(is.empty.bg)) {
    bc@background <- matrix(ncol = 0, nrow = 0)
  } else {
    bg <- lapply(bg, FUN = function(y) y[, common.cells, drop = FALSE])
    background <- as.matrix(do.call("rbind", bg[!is.empty.bg]))
    rownames(background) <- gsub("bc[1|2]\\.", "", rownames(background))
    bc@background <- background[unique(rownames(background)), , drop = FALSE]
  } 
  return(bc)
}


#' @title Creates a new beyondcell object
#' @description This function creates a new \code{\link[beyondcell]{beyondcell}}
#' object.
#' @name CreatebcObject
#' @import scales
#' @param bc \code{beyondcell} object.
#' @details This function creates a new \code{beyondcell} object by using the
#' normalized BCS as the original \code{@@data}. Switch points are
#' recomputed and \code{@@regression} is restarted. The \code{@@expr.matrix} and
#' \code{@@meta.data} slots are subsetted to keep only those cells present in
#' the new \code{@@data} slot.
#' @return A new \code{beyondcell} object.
#' @examples
#' @export

CreatebcObject  <- function(bc) {
  # --- Checks ---
  # Check that bc is a beyondcell object.
  if (class(bc) != "beyondcell") stop('bc must be a beyondcell object.')
  # --- Code ---
  # Make @data equal to @normalized.
  bc@data <- bc@normalized
  # Recompute switch points.
  bc@switch.point <- SwitchPoint(bc)
  # Subset those cells in @expr.matrix that have been removed from @data.
  bc@expr.matrix <- bc@expr.matrix[, colnames(bc@data), drop = FALSE]
  # Subset those cells in @meta.data that have been removed from @data.
  bc@meta.data <- bc@meta.data[colnames(bc@data), , drop = FALSE]
  # Restore @regression.
  bc@regression <- list(order = rep("", 2), vars = NULL,
                        order.background = rep("", 2))
  # Restore @background.
  bc@background <- matrix(ncol = 0, nrow = 0)
  return(bc)
}
