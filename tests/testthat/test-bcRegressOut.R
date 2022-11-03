# PBMC data.
pbmc.data <- Seurat::Read10X("../testdata/single-cell/", gene.column = 1)

# Seurat object.
pbmc.raw <- Seurat::CreateSeuratObject(counts = pbmc.data, project = "pbmc3k",
                                       min.cells = 3, min.features = 20)

# Normalize.
pbmc <- Seurat::NormalizeData(pbmc.raw, normalization.method = "LogNormalize",
                              scale.factor = 10000)

# Geneset objects.
gs100 <- GenerateGenesets("../testdata/gmt/correct100.gmt")
gs10 <- GenerateGenesets("../testdata/gmt/correct10.gmt")
gs20 <- GenerateGenesets("../testdata/gmt/correct20.gmt")

# Beyondcell object.
bc.pbmc <- bcScore(pbmc, gs = gs100, expr.thres = 0.1)
bc.pbmc10 <- bcScore(pbmc, gs = gs10, expr.thres = 0.1)
bc.pbmc20 <- bcScore(pbmc, gs = gs20, expr.thres = 0.1)

# Beyondcell object that has been already regressed.
bc.reg1 <- bc.reg2 <- bc.pbmc
bc.reg1@regression$order <- c("regression", "")
#bc.reg2 <- bcUMAP(bc.reg2, pc = 10, add.DSS = TRUE)
bc.reg2@regression$order <- c("subset", "regression")
bc.reg2@regression$order.background <- c("subset", "regression")

# Beyondcell object that has been already regressed and subsetted.
bc.reg.sub <- bc.pbmc
bc.reg.sub@regression$order <- c("regression", "subset")

# Corrupt beyondcell objects
bc.corrupt1 <- bc.corrupt2 <- bc.corrupt3 <- bc.corrupt4 <- bc.pbmc
bc.corrupt5 <- bc.reg2
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
    bcRegressOut(bc.pbmc, vars.to.regress = "a"),
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
    bcRegressOut(bc.pbmc, vars.to.regress = "nFeature_RNA", k.neighbors = "a"),
    'k.neighbors must be a positive integer.'
  )
  testthat::expect_error(
    bcRegressOut(bc.pbmc, vars.to.regress = "nFeature_RNA", k.neighbors = 1:2),
    'k.neighbors must be a positive integer.'
  )
  testthat::expect_error(
    bcRegressOut(bc.pbmc, vars.to.regress = "nFeature_RNA", k.neighbors = 1.5),
    'k.neighbors must be a positive integer.'
  )
  testthat::expect_error(
    bcRegressOut(bc.pbmc, vars.to.regress = "nFeature_RNA", k.neighbors = -2),
    'k.neighbors must be a positive integer.'
  )
  ### Check add.DSS.
  testthat::expect_error(
    bcRegressOut(bc.pbmc, vars.to.regress = "nFeature_RNA", 
                 add.DSS = rep(TRUE, 2)),
    'add.DSS must be TRUE or FALSE.'
  )
  testthat::expect_error(
    bcRegressOut(bc.pbmc, vars.to.regress = "nFeature_RNA", add.DSS = 1),
    'add.DSS must be TRUE or FALSE.'
  )
  testthat::expect_error(
    bcRegressOut(bc.pbmc10, vars.to.regress = "nFeature_RNA", add.DSS = FALSE),
    paste('Only 10 drug signatures (excluding pathways) are present in the bc', 
          'object, please set add.DSS = TRUE.'), fixed = TRUE
  )
})

# Test warnings.
testthat::test_that("warnings", {
  ### Check regression and subset order.
  testthat::expect_equal(
    testthat::capture_warning(
      bcRegressOut(bc.corrupt1, vars.to.regress = "nFeature_RNA"),
    )$message,
    'Corrupt beyondcell object. Restoring original object before regressing...'
  )
  testthat::expect_equal(
    testthat::capture_warning(
      bcRegressOut(bc.corrupt2, vars.to.regress = "nFeature_RNA"),
    )$message,
    'Corrupt beyondcell object. Restoring original object before regressing...'
  )
  testthat::expect_equal(
    testthat::capture_warning(
      bcRegressOut(bc.corrupt3, vars.to.regress = "nFeature_RNA"),
    )$message,
    'Corrupt beyondcell object. Restoring original object before regressing...'
  )
  testthat::expect_equal(
    testthat::capture_warning(
      bcRegressOut(bc.corrupt4, vars.to.regress = "nFeature_RNA"),
    )$message,
    'Corrupt beyondcell object. Restoring original object before regressing...'
  )
  testthat::expect_equal(
    testthat::capture_warning(
      bcRegressOut(bc.corrupt5, vars.to.regress = "nFeature_RNA"),
    )$message,
    'Corrupt beyondcell object. Restoring original object before regressing...'
  )
  testthat::expect_equal(
    testthat::capture_warning(
      bcRegressOut(bc.reg1, vars.to.regress = "nFeature_RNA"),
    )$message,
    'bc is an already regressed object.'
  )
  testthat::expect_equal(
    testthat::capture_warning(
      bcRegressOut(bc.reg2, vars.to.regress = "nFeature_RNA"),
    )$message,
    'bc is an already regressed object.'
  )
  ### Check add.DSS.
  testthat::expect_equal(
    testthat::capture_warning(
      bcRegressOut(bc.pbmc20, vars.to.regress = "nFeature_RNA", 
                   add.DSS = FALSE),
    )$message,
    paste('Computing an UMAP reduction for 20 drugs. We recommend to set', 
          'add.DSS = TRUE when the number of signatures (excluding pathways)', 
          'is below or equal to 20.')
  )
})