# --- Functions ---
# Function to capitalize genes (mouse symbol-like).
capitalize <- function(x) {
  x <- tolower(x)
  substr(x, 1, 1) <- toupper(substr(x, 1, 1))
  return(x)
}

# --- Code ---
# PBMC data.
pbmc.data <- Seurat::Read10X("../testdata/")

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

# Geneset object.
ssc <- GetCollection(SSc, n.genes = 100)

# Geneset with "mouse" genes.
ssc.mouse <- ssc
ssc.mouse@genelist <- lapply(ssc.mouse@genelist, FUN = function(x) {
  lapply(x, function(y) capitalize(y))})
ssc.mouse@info <- data.frame()

# Small experiment with low expression cells for 2/3 of signatures.
### Small geneset with 3 signatures (only up).
ngenes <- 50
ssc.lowexpr <- GetCollection(SSc, n.genes = ngenes, mode = "up",
                             filters = list(IDs = c("sig-20879", "sig-20880", 
                                                    "sig-20881")),
                             include.pathways = FALSE)
### Unique genes in all 3 signatures.
all.genes <- unique(unlist(ssc.lowexpr@genelist))
### Subset pbmc.data to keep as many rows as unique genes.
pbmc.data.lowexpr <- pbmc.data[1:length(all.genes), ]
### Change rownames.
rownames(pbmc.data.lowexpr) <- all.genes
### Make all expression values = 0...
pbmc.data.lowexpr[, ] <- 10
### Except for >10% of the genes in signatures 2 and 3 (expression = 10).
pbmc.data.lowexpr[ssc.lowexpr@genelist[[2]][["up"]][1:((0.1*ngenes)+1)], ] <- 10
pbmc.data.lowexpr[ssc.lowexpr@genelist[[3]][["up"]][1:((0.1*ngenes)+1)], ] <- 10
### Create Seurat object.
pbmc.raw.lowexpr <- Seurat::CreateSeuratObject(counts = pbmc.data.lowexpr, 
                                               project = "lowexpr",
                                               min.cells = 3, min.features = 2)
### And normalize.
pbmc.lowexpr <- Seurat::NormalizeData(pbmc.raw.lowexpr, scale.factor = 10000, 
                                      normalization.method = "LogNormalize")

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

# Test warnings.
testthat::test_that("warnings", {
  ### Check when sc is an expression matrix.
  testthat::expect_equal(
    testthat::capture_warning(
      bcScore(mtx, gs = ssc)
    )$message,
    paste('Using count matrix as input. Please, check that this matrix',
          'is normalized and unscaled.')
  )
  ### Check when gene names are not in the same format.
  testthat::expect_equal(
    testthat::capture_warning(
      bcScore(pbmc, gs = ssc.mouse)
    )$message,
    paste('gs genes are capitalized and sc genes are in uppercase. Please', 
          'check your Seurat object and translate the genes if necessary.')
  )
  ### Check signatures for which no cells pass the expr.thres.
  testthat::expect_equal(
    testthat::capture_warning(
      bcScore(pbmc.lowexpr, gs = ssc.lowexpr)
    )$message,
    paste('The following signatures have no cells that pass the expr.thres and',
          'will be removed: "sig-20880", "sig-20881".')
  )
})