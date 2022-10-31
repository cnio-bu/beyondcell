# PBMC data.
pbmc.data <- Seurat::Read10X("../testdata/")

# Seurat object.
pbmc.raw <- Seurat::CreateSeuratObject(counts = pbmc.data, project = "pbmc3k",
                                       min.cells = 3, min.features = 20)

# Normalize
pbmc <- Seurat::NormalizeData(pbmc.raw, normalization.method = "LogNormalize",
                              scale.factor = 10000)

# Fake assay.
pbmc.fake <- pbmc
pbmc.fake[["ADT"]] <- Seurat::CreateAssayObject(counts = pbmc.data)
Seurat::DefaultAssay(pbmc.fake) <- "ADT"

# Geneset object.
ssc <- GetCollection(SSc, n.genes = 100)

# Test errors.
testthat::test_that("errors", {
  ### Check sc.
  testthat::expect_error(
    bcScore(1:10, gs = ssc),
    'sc must be either a Seurat object or a single-cell expression matrix.',
    fixed = TRUE
  )
  testthat::expect_error(
    bcScore(pbmc.raw, gs = ssc),
    'Default assay must include a normalized data (@data) slot.',
    fixed = TRUE
  )
  testthat::expect_error(
    bcScore(pbmc.fake, gs = ssc),
    'Seurat default assay must be either RNA, Spatial or SCT.'
  )
  ### Check gs.
  testthat::expect_error(
    bcScore(pbmc, gs = ssc@genelist),
    'gs must be a geneset object.'
  )
  ### Check expr.thres.
  testthat::expect_error(
    bcScore(pbmc, gs = ssc, expr.thres = "a"),
    'expr.thres must be a positive number between 0 and 1.'
  )
  testthat::expect_error(
    bcScore(pbmc, gs = ssc, expr.thres = 1:2),
    'expr.thres must be a positive number between 0 and 1.'
  )
  testthat::expect_error(
    bcScore(pbmc, gs = ssc, expr.thres = 1.5),
    'expr.thres must be a positive number between 0 and 1.'
  )
  testthat::expect_error(
    bcScore(pbmc, gs = ssc, expr.thres = -2),
    'expr.thres must be a positive number between 0 and 1.'
  )
})
