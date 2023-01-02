# --- Functions ---
# Function to subset a geneset object by mode and/or length.
subset.geneset <- function(gs, mode, n) {
  gs@genelist <- gs@genelist[1:n]
  gs@genelist <- lapply(gs@genelist, FUN = function(x) x[mode])
  gs@n.genes <- median(sapply(gs@genelist, FUN = function(y) length(unlist(y))))
  gs@mode <- mode
  return(gs)
}

# Function to compute the scaled BCS.
scale.score <- function(normalized) {
  scaled <- round(t(apply(normalized, 1, scales::rescale, to = 0:1)), 
                  digits = 2)
  return(scaled)
}

# Function to compute the switch point when mode = c("up", "down").
sp.up.down <- function(normalized, scaled) {
  sapply(1:nrow(normalized), FUN = function(i) {
    scores <- normalized[i, ]
    if (all(na.omit(scores) > 0)) return(0)
    else if (all(na.omit(scores) < 0)) return(1)
    else if (any(na.omit(scores) == 0)) {
      idx <- which(scores == 0 & !is.na(scores))[1]
      round(scaled[i, idx], digits = 2)
    } else {
      min.value <- min(normalized[i, scores > 0 & !is.na(scores)])
      min.idx <- which(normalized[i, ] == min.value & !is.na(scaled[i, ]))[1]
      max.value <- max(normalized[i, scores < 0 & !is.na(scores)])
      max.idx <- which(normalized[i, ] == max.value & !is.na(scaled[i, ]))[1]
      round(mean(scaled[i, c(min.idx, max.idx)]), digits = 2)
    }
  })
}

# --- Code ---
# Seed.
set.seed(1)

# PBMC data.
pbmc.data <- Seurat::Read10X("../testdata/single-cell/", gene.column = 1)

# Seurat object.
pbmc.raw <- Seurat::CreateSeuratObject(counts = pbmc.data, project = "pbmc80",
                                       min.cells = 3, min.features = 20)

# Normalize.
pbmc <- Seurat::NormalizeData(pbmc.raw, normalization.method = "LogNormalize",
                              scale.factor = 10000)

# Geneset objects.
gs100 <- GenerateGenesets("../testdata/gmt/correct100.gmt")
gs10 <- GenerateGenesets("../testdata/gmt/correct10.gmt")

gs20 <- subset.geneset(gs100, mode = c("up", "down"), n = 20)
gs100up <- subset.geneset(gs100, mode = "up", n = length(gs100@genelist))
gs100down <- subset.geneset(gs100, mode = "down", n = length(gs100@genelist))

# Beyondcell objects.
bc.object <- suppressWarnings(bcScore(pbmc, gs = gs100, expr.thres = 0.26))
bc.object.bg <- suppressWarnings(
  bcUMAP(bc.object, pc = 2, add.DSS = TRUE)
)
indexes <- length(bc.object.bg@background)
bc.object.bg@background[sample(indexes, size = 100, replace = FALSE)] <- NaN

bc.object.up <- bcScore(pbmc, gs = gs100up, expr.thres = 0.1)
bc.object.down <- bcScore(pbmc, gs = gs100down, expr.thres = 0.1)

bc.object.complete <- bcScore(pbmc, gs = gs100, expr.thres = 0)
bc.object.complete.bg <- suppressWarnings(
  bcUMAP(bc.object.complete, pc = 2, add.DSS = TRUE)
)

bc.object.norm.complete <- bc.object
bc.object.norm.complete.bg <- bc.object.bg
bc.object.norm.complete@normalized <- 
  bc.object.norm.complete.bg@normalized <- bc.object.complete@normalized

bc.object10 <- bcScore(pbmc, gs = gs10, expr.thres = 0.1)
bc.object20 <- bcScore(pbmc, gs = gs20, expr.thres = 0.1)

# Beyondcell object that has already been subsetted.
bc.sub <- bc.object
bc.sub@regression$order <- c("subset", "")

# Beyondcell object with background matrix that has already been regressed.
bc.reg.bg <- bc.object.bg
bc.reg.bg@regression$order <- c("regression", "")
bc.reg.bg@regression$order.background <- c("regression", "")

# Beyondcell object that has already been subsetted and regressed.
bc.sub.reg <- bc.object
bc.sub.reg@regression$order <- c("subset", "regression")

# Beyondcell object that has already been regressed and subsetted.
bc.reg.sub <- bc.object
bc.reg.sub@regression$order <- c("regression", "subset")

# Beyondcell object with background matrix that has been already regressed and 
# subsetted.
bc.reg.sub.bg <- bc.object.bg
bc.reg.sub.bg@regression$order <- c("regression", "subset")
bc.reg.sub.bg@regression$order.background <- c("regression", "subset")
bc.reg.sub.bg@background <- bc.reg.sub.bg@background[1:10, ]

# Corrupt beyondcell objects.
bc.corrupt1 <- bc.corrupt2 <- bc.corrupt3 <- bc.corrupt4 <- bc.object
bc.corrupt5 <- bc.reg.bg
bc.corrupt1@regression$order <- c("a", "")
bc.corrupt2@regression$order <- c("", "subset")
bc.corrupt3@regression$order <- rep("", 3)
bc.corrupt4@regression$order <- rep("subset", 2)
bc.corrupt5@regression$order <- c("subset", "")

# Beyondcell object with only one signature.
gs1 <- GetCollection(SSc, n.genes = 100, mode = c("up", "down"),
                     filters = list(IDs = SSc@info$IDs[1]),
                     include.pathways = FALSE)
bc1 <- bcScore(pbmc, gs = gs1, expr.thres = 0.1)

# Test errors.
testthat::test_that("errors", {
  ### Check bc.
  testthat::expect_error(
    bcRegressOut(letters[1:5], vars.to.regress = "a"),
    'bc must be a beyondcell object.'
  )
  ### Check vars.to.regress.
  testthat::expect_error(
    bcRegressOut(bc.object, vars.to.regress = "a"),
    'vars.to.regress not found.'
  )
  ### Check regression and subset order.
  testthat::expect_error(
    bcRegressOut(bc.reg.sub, vars.to.regress = "nFeature_RNA"),
    paste('bc was previously regressed and then subsetted. Please run',
          'bcSubset() on a new beyondcell object created with', 
          'CreatebcObject(bc).'), fixed = TRUE
  )
  ### Check add.DSS.
  testthat::expect_error(
    bcRegressOut(bc.object, vars.to.regress = "nFeature_RNA", 
                 add.DSS = rep(TRUE, 2)),
    'add.DSS must be TRUE or FALSE.'
  )
  testthat::expect_error(
    bcRegressOut(bc.object, vars.to.regress = "nFeature_RNA", add.DSS = 1),
    'add.DSS must be TRUE or FALSE.'
  )
  testthat::expect_error(
    bcRegressOut(bc.object10, vars.to.regress = "nFeature_RNA", add.DSS = FALSE),
    paste('Only 10 drug signatures (excluding pathways) are present in the bc', 
          'object, please set add.DSS = TRUE.'), fixed = TRUE
  )
  ### Check k.neighbors.
  n.complete.norm <- sum(complete.cases(t(bc.object@normalized)))
  n.complete.bg <- sum(complete.cases(t(bc.object.bg@background)))
  bc.object.norm.complete.thres1 <- bc.object.norm.complete
  bc.object.norm.complete.thres1@thres <- 1
  testthat::expect_error(
    bcRegressOut(bc.object, vars.to.regress = "nFeature_RNA", k.neighbors = "a"),
    'k.neighbors must be numeric.'
  )
  testthat::expect_error(
    bcRegressOut(bc.object, vars.to.regress = "nFeature_RNA", k.neighbors = 1:2),
    'k.neighbors must be a positive integer.'
  )
  testthat::expect_error(
    bcRegressOut(bc.object, vars.to.regress = "nFeature_RNA", k.neighbors = 1.5),
    'k.neighbors must be a positive integer.'
  )
  testthat::expect_error(
    bcRegressOut(bc.object, vars.to.regress = "nFeature_RNA", k.neighbors = 0),
    'k.neighbors must be a positive integer.'
  )
  testthat::expect_error(
    bcRegressOut(bc.object, vars.to.regress = "nFeature_RNA", k.neighbors = -2),
    'k.neighbors must be a positive integer.'
  )
  testthat::expect_error(
    bcRegressOut(bc.object, vars.to.regress = "nFeature_RNA", 
                 k.neighbors = n.complete.norm, add.DSS = FALSE),
    paste0('k.neighbors must be lower than the number of complete cases in ', 
           '@normalized slot: ', n.complete.norm, '.')
  )
  testthat::expect_error(
    bcRegressOut(bc.object, vars.to.regress = "nFeature_RNA", 
                 k.neighbors = n.complete.norm + 1, add.DSS = FALSE),
    paste0('k.neighbors must be lower than the number of complete cases in ', 
           '@normalized slot: ', n.complete.norm, '.')
  )
  testthat::expect_error(
    bcRegressOut(bc.object.norm.complete.bg, vars.to.regress = "nFeature_RNA", 
                 k.neighbors = n.complete.bg, add.DSS = FALSE),
    paste0('k.neighbors must be lower than the number of complete cases in ', 
           '@background slot: ', n.complete.bg, '.')
  )
  testthat::expect_error(
    bcRegressOut(bc.object.norm.complete.thres1, 
                 vars.to.regress = "nFeature_RNA", 
                 k.neighbors = n.complete.norm, add.DSS = TRUE),
    paste0('k.neighbors must be lower than the total number of complete ',
           'cases in @normalized and @background slots: ', n.complete.norm, '.')
  )
  ### Check that bcRegressOut does not throw any error with only one signature.
  testthat::expect_no_error(
    bcRegressOut(bc1, vars.to.regress = "nFeature_RNA", add.DSS = TRUE)
  )
})

# Test warnings.
testthat::test_that("warnings", {
  ### Check regression and subset order.
  testthat::expect_warning(
    bcRegressOut(bc.corrupt1, vars.to.regress = "nFeature_RNA"),
    'Corrupt beyondcell object. Restoring original object before regressing...'
  )
  testthat::expect_warning(
    bcRegressOut(bc.corrupt2, vars.to.regress = "nFeature_RNA"),
    'Corrupt beyondcell object. Restoring original object before regressing...'
  )
  testthat::expect_warning(
    bcRegressOut(bc.corrupt3, vars.to.regress = "nFeature_RNA"),
    'Corrupt beyondcell object. Restoring original object before regressing...'
  )
  testthat::expect_warning(
    bcRegressOut(bc.corrupt4, vars.to.regress = "nFeature_RNA"),
    'Corrupt beyondcell object. Restoring original object before regressing...'
  )
  testthat::expect_warning(
    bcRegressOut(bc.corrupt5, vars.to.regress = "nFeature_RNA"),
    'Corrupt beyondcell object. Restoring original object before regressing...'
  )
  testthat::expect_warning(
    bcRegressOut(bc.sub.reg, vars.to.regress = "nFeature_RNA"),
    'bc is an already regressed object.'
  )
  testthat::expect_warning(
    bcRegressOut(bc.reg.bg, vars.to.regress = "nFeature_RNA"),
    'bc is an already regressed object.'
  )
  ### Check add.DSS.
  testthat::expect_warning(
    bcRegressOut(bc.object20, vars.to.regress = "nFeature_RNA", add.DSS = FALSE),
    paste('Computing an UMAP reduction for 20 drugs. We recommend to set', 
          'add.DSS = TRUE when the number of signatures (excluding pathways)', 
          'is below or equal to 20.'), fixed = TRUE
  )
})

# Test messages.
testthat::test_that("messages", {
  ### Check the printed messages when using a beyondcell object without
  ### background matrix and add.DSS = FALSE as inputs.
  testthat::expect_equal(
    testthat::capture_messages(
      bcRegressOut(bc.object, vars.to.regress = "nFeature_RNA", add.DSS = FALSE)
    ),
    c(paste('DSS background not computed. The imputation will be computed with',
            'just the drugs (not pathways) in the beyondcell object.\n'),
      'Imputing normalized BCS...\n',
      'Regressing scores...\n')
  )
  ### Check the printed messages when using a beyondcell object without
  ### background matrix and add.DSS = TRUE as inputs.
  testthat::expect_equal(
    testthat::capture_messages(
      bcRegressOut(bc.object, vars.to.regress = "nFeature_RNA", add.DSS = TRUE)
    ),
    c('Computing background BCS using DSS signatures...\n',
      'Imputing normalized BCS...\n',
      'Regressing scores...\n',
      'Imputing background BCS...\n',
      'Regressing background BCS...\n')
  )
  testthat::expect_equal(
    testthat::capture_messages(
      suppressWarnings(
        bcRegressOut(bc.sub.reg, vars.to.regress = "nFeature_RNA", 
                     add.DSS = TRUE)
      )
    ),
    c('Computing background BCS using DSS signatures...\n',
      'Imputing normalized BCS...\n',
      'Regressing scores...\n',
      'Imputing background BCS...\n',
      'Regressing background BCS...\n')
  )
  ### Check the printed messages when using a beyondcell object with
  ### background matrix and add.DSS = FALSE as inputs.
  testthat::expect_equal(
    testthat::capture_messages(
      bcRegressOut(bc.object.bg, vars.to.regress = "nFeature_RNA", 
                   add.DSS = FALSE)
    ),
    c(paste('DSS background not computed. The imputation will be computed with',
            'just the drugs (not pathways) in the beyondcell object.\n'),
      'Imputing normalized BCS...\n',
      'Regressing scores...\n',
      'Removing @reductions slot...\n',
      'Removing therapeutic clusters...\n',
      'Imputing background BCS...\n',
      'Regressing background BCS...\n')
  )
  ### Check the printed messages when using a beyondcell object with
  ### background matrix and add.DSS = TRUE as inputs.
  testthat::expect_equal(
    testthat::capture_messages(
      bcRegressOut(bc.object.bg, vars.to.regress = "nFeature_RNA", 
                   add.DSS = TRUE)
    ),
    c('Background BCS already computed. Skipping this step.\n',
      'Imputing normalized BCS...\n',
      'Regressing scores...\n',
      'Removing @reductions slot...\n',
      'Removing therapeutic clusters...\n',
      'Imputing background BCS...\n',
      'Regressing background BCS...\n')
  )
  testthat::expect_equal(
    testthat::capture_messages(
      suppressWarnings(
        bcRegressOut(bc.reg.bg, vars.to.regress = "nFeature_RNA", 
                     add.DSS = TRUE)
      )
    ),
    c('Restoring pre-regressed background matrix...\n',
      'Background BCS already computed. Skipping this step.\n',
      'Imputing normalized BCS...\n',
      'Regressing scores...\n',
      'Imputing background BCS...\n',
      'Regressing background BCS...\n')
  )
  ### Check the printed messages when using a complete beyondcell object without
  ### background matrix and add.DSS = FALSE as inputs.
  testthat::expect_equal(
    testthat::capture_messages(
      bcRegressOut(bc.object.complete, vars.to.regress = "nFeature_RNA", 
                   add.DSS = FALSE)
    ),
    c(paste('DSS background not computed. The imputation will be computed with',
            'just the drugs (not pathways) in the beyondcell object.\n'),
      'No NaN values were found in bc@normalized. No imputation needed.\n',
      'Regressing scores...\n')
  )
  ### Check the printed messages when using a complete beyondcell object without
  ### background matrix and add.DSS = TRUE as inputs.
  testthat::expect_equal(
    testthat::capture_messages(
      bcRegressOut(bc.object.complete, vars.to.regress = "nFeature_RNA", 
                   add.DSS = TRUE)
    ),
    c('Computing background BCS using DSS signatures...\n',
      'No NaN values were found in bc@normalized. No imputation needed.\n',
      'Regressing scores...\n',
      'No NaN values were found in bc@background. No imputation needed.\n',
      'Regressing background BCS...\n')
  )
  testthat::expect_equal(
    testthat::capture_messages(
      bcRegressOut(bc.object.norm.complete, vars.to.regress = "nFeature_RNA", 
                   add.DSS = TRUE)
    ),
    c('Computing background BCS using DSS signatures...\n',
      'No NaN values were found in bc@normalized. No imputation needed.\n',
      'Regressing scores...\n',
      'Imputing background BCS...\n',
      'Regressing background BCS...\n')
  )
  ### Check the printed messages when using a complete beyondcell object with
  ### background matrix and add.DSS = FALSE as inputs.
  testthat::expect_equal(
    testthat::capture_messages(
      bcRegressOut(bc.object.complete.bg, add.DSS = FALSE,
                   vars.to.regress = "nFeature_RNA")
    ),
    c(paste('DSS background not computed. The imputation will be computed with',
            'just the drugs (not pathways) in the beyondcell object.\n'),
      'No NaN values were found in bc@normalized. No imputation needed.\n',
      'Regressing scores...\n',
      'Removing @reductions slot...\n',
      'Removing therapeutic clusters...\n',
      'No NaN values were found in bc@background. No imputation needed.\n',
      'Regressing background BCS...\n')
  )
  testthat::expect_equal(
    testthat::capture_messages(
      bcRegressOut(bc.object.norm.complete.bg, add.DSS = FALSE,
                   vars.to.regress = "nFeature_RNA")
    ),
    c(paste('DSS background not computed. The imputation will be computed with',
            'just the drugs (not pathways) in the beyondcell object.\n'),
      'No NaN values were found in bc@normalized. No imputation needed.\n',
      'Regressing scores...\n',
      'Removing @reductions slot...\n',
      'Removing therapeutic clusters...\n',
      'Imputing background BCS...\n',
      'Regressing background BCS...\n')
  )
  ### Check the printed messages when using a complete beyondcell object with
  ### background matrix and add.DSS = TRUE as inputs.
  testthat::expect_equal(
    testthat::capture_messages(
      bcRegressOut(bc.object.complete.bg, add.DSS = TRUE,
                   vars.to.regress = "nFeature_RNA")
    ),
    c('Background BCS already computed. Skipping this step.\n',
      'No NaN values were found in bc@normalized. No imputation needed.\n',
      'Regressing scores...\n',
      'Removing @reductions slot...\n',
      'Removing therapeutic clusters...\n',
      'No NaN values were found in bc@background. No imputation needed.\n',
      'Regressing background BCS...\n')
  )
  testthat::expect_equal(
    testthat::capture_messages(
      bcRegressOut(bc.object.norm.complete.bg, add.DSS = FALSE,
                   vars.to.regress = "nFeature_RNA")
    ),
    c('Background BCS already computed. Skipping this step.\n',
      'No NaN values were found in bc@normalized. No imputation needed.\n',
      'Regressing scores...\n',
      'Removing @reductions slot...\n',
      'Removing therapeutic clusters...\n',
      'Imputing background BCS...\n',
      'Regressing background BCS...\n')
  )
})

# Test values.
testthat::test_that("default values", {
  bc.regressed <- bcRegressOut(bc.object, add.DSS = FALSE,
                               vars.to.regress = "nFeature_RNA")
  bc.regressed.bg <- bcRegressOut(bc.object, add.DSS = TRUE,
                                  vars.to.regress = "nFeature_RNA")
  bc.regressed.up <- bcRegressOut(bc.object.up, add.DSS = FALSE,
                                  vars.to.regress = "nFeature_RNA")
  bc.regressed.down <- bcRegressOut(bc.object.down, add.DSS = FALSE,
                                    vars.to.regress = "nFeature_RNA")
  ### Check that bcRegressOut output is a beyondcell object.
  testthat::expect_s4_class(
    bc.regressed,
    "beyondcell"
  )
  ### Check that the slot @normalized is equal to the slot @data.
  testthat::expect_equal(
    bc.regressed@normalized,
    bc.regressed@data
  )
  ### Check that the slot @scaled is a matrix.
  testthat::expect_equal(
    class(bc.regressed@scaled),
    c("matrix", "array")
  )
  ### Check that the slot @normalized is a matrix.
  testthat::expect_equal(
    class(bc.regressed@normalized),
    c("matrix", "array")
  )
  ### Check that the slot @background is a matrix.
  testthat::expect_equal(
    class(bc.regressed.bg@background),
    c("matrix", "array")
  )
  ### Check that the matrices in slots @scaled, @normalized and @data have the 
  ### same dimensions.
  testthat::expect_true(
    identical(dim(bc.regressed@scaled), dim(bc.regressed@normalized))
  )
  ### Check that the @background matrix has the same number of columns as the 
  ### other matrices.
  testthat::expect_true(
    identical(ncol(bc.regressed.bg@scaled), ncol(bc.regressed.bg@background))
  )
  ### Check that the matrices in slots @scaled, @normalized and @data have the 
  ### same dimnames.
  testthat::expect_true(
    identical(dimnames(bc.regressed@scaled), dimnames(bc.regressed@normalized))
  )
  ### Check that the @background matrix has the same colnames as the other 
  ### matrices.
  testthat::expect_true(
    identical(colnames(bc.regressed.bg@scaled), 
              colnames(bc.regressed.bg@background))
  )
  ### Check that the matrices do not contain NaN values.
  testthat::expect_false(
    any(is.na(bc.regressed@scaled))
  )
  testthat::expect_false(
    any(is.na(bc.regressed@normalized))
  )
  testthat::expect_false(
    any(is.na(bc.regressed.bg@background))
  )
  ### Check that the values in the slot @scaled range between 0 and 1.
  testthat::expect_equal(
    c(min(bc.regressed@scaled, na.rm = TRUE), 
      max(bc.regressed@scaled, na.rm = TRUE)),
    0:1
  )
  ### Check that the values in the slot @scaled are different between input and
  ### output.
  testthat::expect_false(
    isTRUE(all.equal(bc.regressed@scaled, bc.object@scaled))
  )
  ### Check the values of the slot @scaled.
  testthat::expect_equal(
    bc.regressed@scaled, scale.score(bc.regressed@normalized)
  )
  testthat::expect_equal(
    bc.regressed.up@scaled, scale.score(bc.regressed.up@normalized)
  )
  testthat::expect_equal(
    bc.regressed.down@scaled, scale.score(bc.regressed.down@normalized)
  )
  ### Check that the values in the slot @normalized range between -Inf and +Inf.
  testthat::expect_true(
    min(bc.regressed@normalized, na.rm = TRUE) < 0
  )
  testthat::expect_true(
    max(bc.regressed@normalized, na.rm = TRUE) > 0
  )
  ### Check that the values in the slot @normalized are different between input 
  ### and output.
  testthat::expect_false(
    isTRUE(all.equal(bc.regressed@normalized, bc.object@normalized))
  )
  ### Check that the slot @switch.point is a numeric vector.
  testthat::expect_equal(
    class(bc.regressed@switch.point),
    "numeric"
  )
  ### Check that the @switch.point vector doesn't contain NA values.
  testthat::expect_true(
    all(!is.na(bc.regressed@switch.point))
  )
  ### Check that the values in the slot @switch.point range between 0 and 1.
  testthat::expect_true(
    min(bc.regressed@switch.point) >= 0
  )
  testthat::expect_true(
    max(bc.regressed@switch.point) <= 1
  )
  ### Check that the values in the slot @switch.point are different between 
  ### input and output.
  testthat::expect_false(
    isTRUE(all.equal(bc.regressed@switch.point, bc.object@switch.point))
  )
  ### Check the values of the slot @switch.point.
  testthat::expect_equal(
    unname(bc.regressed@switch.point),
    sp.up.down(bc.regressed@normalized, bc.regressed@scaled)
  )
  testthat::expect_equal(
    unique(bc.regressed.up@switch.point),
    0
  )
  testthat::expect_equal(
    unique(bc.regressed.down@switch.point),
    1
  )
  ### Check that the slot @ranks is empty.
  testthat::expect_equal(
    bc.regressed@ranks,
    list()
  )
  ### Check that the slot @expr.matrix contains the input matrix.
  testthat::expect_equal(
    bc.regressed@expr.matrix,
    bc.object@expr.matrix
  )
  # Check that the slot @meta.data contains the input metadata except for TCs 
  # column(s).
  TCs <- which(startsWith(colnames(bc.object.bg@meta.data), "bc_clusters_res."))
  testthat::expect_false(
    any(startsWith(colnames(bc.regressed@meta.data), "bc_clusters_res."))
  )
  testthat::expect_false(
    any(startsWith(colnames(bc.regressed.bg@meta.data), "bc_clusters_res."))
  )
  testthat::expect_equal(
    bc.regressed@meta.data,
    bc.object@meta.data
  )
  testthat::expect_equal(
    bc.regressed.bg@meta.data,
    bc.object.bg@meta.data[, -c(TCs)]
  )
  ### Check the values of the slot @SeuratInfo.
  testthat::expect_equal(
    bc.regressed@SeuratInfo,
    bc.object@SeuratInfo
  )
  ### Check the values of slot @background.
  testthat::expect_equal(
    bc.regressed@background,
    bc.object@background
  )
  ### Check that the slot @reductions is empty.
  testthat::expect_equal(
    bc.regressed@reductions,
    list()
  )
  testthat::expect_equal(
    bc.regressed.bg@reductions,
    list()
  )
  ### Check the values of slot @regression.
  ordering <- c("order", "vars", "order.background")
  testthat::expect_equal(
    bc.regressed@regression[ordering],
    list(order = c("regression", ""), vars = "nFeature_RNA", 
         order.background = rep("", 2))
  )
  testthat::expect_equal(
    bc.regressed.bg@regression[ordering],
    list(order = c("regression", ""), vars = "nFeature_RNA", 
         order.background = c("regression", ""))
  )
  testthat::expect_equal(
    bcRegressOut(bc.sub, add.DSS = TRUE,
                 vars.to.regress = "nFeature_RNA")@regression[ordering],
    list(order = c("subset", "regression"), vars = "nFeature_RNA", 
         order.background = c("subset", "regression"))
  )
  testthat::expect_equal(
    suppressWarnings(
      bcRegressOut(bc.reg.bg, add.DSS = FALSE,
                   vars.to.regress = "nFeature_RNA")@regression[ordering]
    ),
    list(order = c("regression", ""), vars = "nFeature_RNA", 
         order.background = c("regression", ""))
  )
  testthat::expect_equal(
    suppressWarnings(
      bcRegressOut(bc.sub.reg, add.DSS = FALSE, 
                   vars.to.regress = "nFeature_RNA")@regression[ordering]
    ),
    list(order = c("subset", "regression"), vars = "nFeature_RNA", 
         order.background = rep("", 2))
  )
  testthat::expect_equal(
    suppressWarnings(
      bcRegressOut(bc.corrupt1, add.DSS = TRUE,
                   vars.to.regress = "nFeature_RNA")@regression[ordering]
    ),
    list(order = c("regression", ""), vars = "nFeature_RNA", 
         order.background = c("regression", ""))
  )
  ### Check that the slot @n.genes is equal to the slot @n.genes in the input 
  ### beyondcell.
  testthat::expect_equal(
    bc.regressed@n.genes,
    bc.object@n.genes
  )
  ### Check that the slot @mode is equal to the slot @mode in the input 
  ### beyondcell.
  testthat::expect_equal(
    bc.regressed@mode,
    bc.object@mode
  )
  testthat::expect_equal(
    bc.regressed.up@mode,
    bc.object.up@mode
  )
  testthat::expect_equal(
    bc.regressed.down@mode,
    bc.object.down@mode
  )
  ### Check that the slot @mode is equal to the slot @mode in the input 
  ### beyondcell.
  testthat::expect_equal(
    unique(c(bc.regressed@thres, bc.regressed.up@thres, 
             bc.regressed.down@thres)),
    unique(c(bc.object@thres, bc.object.up@thres, bc.object.down@thres))
  )
})