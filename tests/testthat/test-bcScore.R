# --- Functions ---
# Function to capitalize genes (mouse symbol-like).
capitalize <- function(x) {
  x <- tolower(x)
  substr(x, 1, 1) <- toupper(substr(x, 1, 1))
  return(x)
}

# --- Code ---
# PBMC data.
pbmc.data <- Seurat::Read10X("../testdata/single-cell/", gene.column = 1)

# Seurat object.
pbmc.raw <- Seurat::CreateSeuratObject(counts = pbmc.data, project = "pbmc3k",
                                       min.cells = 3, min.features = 20)

# Normalize.
pbmc <- Seurat::NormalizeData(pbmc.raw, normalization.method = "LogNormalize",
                              scale.factor = 10000)

# Fake assay.
pbmc.fake <- pbmc
pbmc.fake[["ADT"]] <- Seurat::CreateAssayObject(counts = pbmc.data)
Seurat::DefaultAssay(pbmc.fake) <- "ADT"

# Single-cell expression matrix.
mtx <- as.matrix(pbmc@assays$RNA@data)

# Geneset objects.
gs10 <- GenerateGenesets("../testdata/gmt/correct10.gmt")
gs.warning <- GenerateGenesets("../testdata/gmt/score_warning10.gmt")

# Geneset with "mouse" genes.
gs.mouse <- gs10
gs.mouse@genelist <- lapply(gs.mouse@genelist, FUN = function(x) {
  lapply(x, function(y) capitalize(y))
})
gs.mouse@info <- data.frame()

# Test errors.
testthat::test_that("errors", {
  ### Check sc.
  testthat::expect_error(
    bcScore(1:10, gs = gs10),
    'sc must be either a Seurat object or a single-cell expression matrix.'
  )
  testthat::expect_error(
    bcScore(pbmc.raw, gs = gs10),
    'Default assay must include a normalized data (@data) slot.',
    fixed = TRUE
  )
  testthat::expect_error(
    bcScore(pbmc.fake, gs = gs10),
    'Seurat default assay must be either RNA, Spatial or SCT.'
  )
  ### Check gs.
  testthat::expect_error(
    bcScore(pbmc, gs = gs10@genelist),
    'gs must be a geneset object.'
  )
  ### Check expr.thres.
  testthat::expect_error(
    bcScore(pbmc, gs = gs10, expr.thres = "a"),
    'expr.thres must be a positive number between 0 and 1.'
  )
  testthat::expect_error(
    bcScore(pbmc, gs = gs10, expr.thres = 1:2),
    'expr.thres must be a positive number between 0 and 1.'
  )
  testthat::expect_error(
    bcScore(pbmc, gs = gs10, expr.thres = 1.5),
    'expr.thres must be a positive number between 0 and 1.'
  )
  testthat::expect_error(
    bcScore(pbmc, gs = gs10, expr.thres = -2),
    'expr.thres must be a positive number between 0 and 1.'
  )
})

# Test warnings.
testthat::test_that("warnings", {
  ### Check when sc is an expression matrix.
  testthat::expect_warning(bcScore(mtx, gs = gs10),
    paste('Using count matrix as input. Please, check that this matrix',
          'is normalized and unscaled.')
  )
  ### Check when gene names are not in the same format.
  testthat::expect_warning(
    bcScore(pbmc, gs = gs.mouse),
    paste('gs genes are capitalized and sc genes are in uppercase. Please', 
          'check your Seurat object and translate the genes if necessary.')
  )
  ### Check signatures for which no cells pass the expr.thres.
  testthat::expect_warning(
    bcScore(pbmc, gs = gs.warning),
    paste('The following signatures have no cells that pass the expr.thres and',
          'will be removed: sig-20965.')
  )
})