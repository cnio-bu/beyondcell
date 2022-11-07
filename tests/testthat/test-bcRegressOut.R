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
gs20 <- GenerateGenesets("../testdata/gmt/correct20.gmt")

# Beyondcell objects.
bc.object <- bcScore(pbmc, gs = gs100, expr.thres = 0.1)
#bc.object.bg <- bcUMAP(bc.object, pc = 10, add.DSS = TRUE)
bc.object.bg <- bc.object # Delete

bc.object.complete <- bcScore(pbmc, gs = gs100, expr.thres = 0)
#bc.object.complete.bg <- bcUMAP(bc.object.complete, pc = 10, 
#                                        add.DSS = TRUE)
bc.object.complete.bg <- bc.object.complete  # Delete
bc.object10 <- bcScore(pbmc, gs = gs10, expr.thres = 0.1)
bc.object20 <- bcScore(pbmc, gs = gs20, expr.thres = 0.1)

# Beyondcell object with background matrix that has already been regressed.
bc.reg.bg <- bc.object.bg
bc.reg.bg@regression$order <- c("regression", "")
bc.reg.bg@regression$order.bg <- c("regression", "")

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
  ### Check k.neighbors.
  testthat::expect_error(
    bcRegressOut(bc.object, vars.to.regress = "nFeature_RNA", k.neighbors = "a"),
    'k.neighbors must be a positive integer.'
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
    bcRegressOut(bc.object, vars.to.regress = "nFeature_RNA", k.neighbors = -2),
    'k.neighbors must be a positive integer.'
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
          'is below or equal to 20.')
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
      bcRegressOut(bc.sub.reg, vars.to.regress = "nFeature_RNA", add.DSS = TRUE)
    ),
    c('Computing background BCS using DSS signatures...\n',
      'Regressing background BCS...\n',
      'Imputing normalized BCS...\n',
      'Regressing scores...\n',
      'Imputing background BCS...\n',
      'Regressing background BCS...\n')
  )
  testthat::expect_equal(
    testthat::capture_messages(
      bcRegressOut(bc.reg.sub, vars.to.regress = "nFeature_RNA", add.DSS = TRUE)
    ),
    c('Computing background BCS using DSS signatures...\n',
      'Regressing background BCS...\n',
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
      'Imputing background BCS...\n',
      'Regressing background BCS...\n')
  )
  testthat::expect_equal(
    testthat::capture_messages(
      bcRegressOut(bc.reg.sub.bg, vars.to.regress = "nFeature_RNA", 
                   add.DSS = TRUE)
    ),
    c('Computing background BCS using DSS signatures...\n',
      'Regressing background BCS...\n',
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
      'Regressing scores...\n')
  )
  ### Check the printed messages when using a complete beyondcell object with
  ### background matrix and add.DSS = FALSE as inputs.
  testthat::expect_equal(
    testthat::capture_messages(
      bcRegressOut(bc.object.complete.bg, add.DSS = TRUE,
                   vars.to.regress = "nFeature_RNA")
    ),
    c('Background BCS already computed. Skipping this step.\n',
      'No NaN values were found in bc@normalized. No imputation needed.\n',
      'Regressing scores...\n',
      'No NaN values were found in bc@background. No imputation needed.\n',
      'Regressing background BCS...\n')
  )
})