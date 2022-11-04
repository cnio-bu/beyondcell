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
pbmc.raw <- Seurat::CreateSeuratObject(counts = pbmc.data, project = "pbmc80",
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
gs10up <- GenerateGenesets("../testdata/gmt/correct10up.gmt")
gs10down <- GenerateGenesets("../testdata/gmt/correct10down.gmt")
gs.warning <- GenerateGenesets("../testdata/gmt/score_warning10.gmt")

# Geneset with "mouse" genes.
gs.mouse <- gs10
gs.mouse@genelist <- lapply(gs.mouse@genelist, FUN = function(x) {
  lapply(x, function(y) capitalize(y))
})
gs.mouse@info <- data.frame()

# Visium data with SCT normalization.
image <- Seurat::Read10X_Image(image.dir = "../testdata/visium/", 
                               image.name = "tissue_lowres_image.png")
pbmc.visium <- Seurat::Load10X_Spatial(data.dir = "../testdata/visium/", 
                                       filename = "matrix.h5", image = image)
pbmc.visium <- Seurat::SCTransform(pbmc.visium, assay = "Spatial", 
                                   verbose = FALSE, new.assay.name = "pbmc80")

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

# Test messages.
testthat::test_that("messages", {
  ### Check initial message.
  Seurat::DefaultAssay(pbmc) <- "RNA"
  testthat::expect_message(
    bcScore(pbmc, gs = gs10),
    'Using RNA assay as input.'
  )
  Seurat::DefaultAssay(pbmc.visium) <- "SCT"
  testthat::expect_message(
    bcScore(pbmc.visium, gs = gs10),
    'Using SCT assay as input.'
  )
  ### Check final message.
  testthat::expect_message(
    bcScore(pbmc, gs = gs10, expr.thres = 0),
    'There are 80/80 cells without missing values in your beyondcell object.'
  )
})

# Test values.
testthat::test_that("default values", {
  bc.object <- bcScore(pbmc, gs = gs10, expr.thres = 0.1)
  bc.object.up <- bcScore(pbmc, gs = gs10up, expr.thres = 0.1)
  bc.object.down <- bcScore(pbmc, gs = gs10down, expr.thres = 0.1)
  ### Test that bcScore output is a beyondcell object.
  testthat::expect_s4_class(
    bc.object,
    "beyondcell"
  )
  ### Check the values of the slot @mode.
  testthat::expect_equal(
    bc.object@mode,
    c("up", "down")
  )
  testthat::expect_equal(
    bc.object.up@mode,
    "up"
  )
  testthat::expect_equal(
    bc.object.down@mode,
    "down"
  )
  ### Check that the slot @scaled is a matrix.
  testthat::expect_equal(
    class(bc.object@scaled),
    c("matrix", "array")
  )
  ### Check that the slot @normalized is a matrix.
  testthat::expect_equal(
    class(bc.object@normalized),
    c("matrix", "array")
  )
  ### Check that the slot @data is a matrix.
  testthat::expect_equal(
    class(bc.object@data),
    c("matrix", "array")
  )
  ### Check that the matrices in slots @scaled, @normalized and @data have the 
  ### same dimensions.
  testthat::expect_true(
    identical(dim(bc.object@scaled), dim(bc.object@normalized))
  )
  testthat::expect_true(
    identical(dim(bc.object@scaled), dim(bc.object@data))
  )
  ### Check that the matrices in slots @scaled, @normalized and @data have the 
  ### same dimnames.
  testthat::expect_true(
    identical(dimnames(bc.object@scaled), dimnames(bc.object@normalized))
  )
  testthat::expect_true(
    identical(dimnames(bc.object@scaled), dimnames(bc.object@data))
  )
  ### Check that the matrices contain NaN values (not NA).
  testthat::expect_true(
    all(is.nan(bc.object@scaled[is.na(bc.object@scaled)]))
  )
  testthat::expect_true(
    all(is.nan(bc.object@normalized[is.na(bc.object@normalized)]))
  )
  testthat::expect_true(
    all(is.nan(bc.object@data[is.na(bc.object@data)]))
  )
  ### Check that the NaN values are at the same indexes.
  testthat::expect_true(
    identical(which(is.nan(bc.object@scaled)), 
              which(is.nan(bc.object@normalized)))
  )
  testthat::expect_true(
    identical(which(is.nan(bc.object@scaled)), which(is.nan(bc.object@data)))
  )
  ### Check that the non NaN values in the slot @scaled range between 0 and 1.
  testthat::expect_equal(
    c(min(bc.object@scaled, na.rm = TRUE), max(bc.object@scaled, na.rm = TRUE)),
    0:1
  )
  ### Check that the non NaN values in the slot @normalized range between -Inf 
  ### and +Inf when mode is c("up", "down").
  testthat::expect_true(
    min(bc.object@normalized, na.rm = TRUE) < 0
  )
  testthat::expect_true(
    max(bc.object@normalized, na.rm = TRUE) > 0
  )
  ### Check that the non NaN values in the slot @normalized range between 0 
  ### and +Inf when mode is "up".
  testthat::expect_true(
    min(bc.object.up@normalized, na.rm = TRUE) >= 0
  )
  testthat::expect_true(
    max(bc.object.up@normalized, na.rm = TRUE) > 0
  )
  ### Check that the non NaN values in the slot @normalized range between -Inf 
  ### and 0 when mode is "down".
  testthat::expect_true(
    min(bc.object.down@normalized, na.rm = TRUE) < 0
  )
  testthat::expect_true(
    max(bc.object.down@normalized, na.rm = TRUE) <= 0
  )
  ### Check that the non NaN values in the slot @data range between -Inf and 
  ### +Inf when mode is c("up", "down").
  testthat::expect_true(
    min(bc.object@data, na.rm = TRUE) < 0
  )
  testthat::expect_true(
    max(bc.object@data, na.rm = TRUE) > 0
  )
  ### Check that the non NaN values in the slot @data range between 0 and 
  ### +Inf when mode is "up".
  testthat::expect_true(
    min(bc.object.up@data, na.rm = TRUE) >= 0
  )
  testthat::expect_true(
    max(bc.object.up@data, na.rm = TRUE) > 0
  )
  ### Check that the non NaN values in the slot @data range between -Inf and 
  ### 0 when mode is "down".
  testthat::expect_true(
    min(bc.object.down@data, na.rm = TRUE) < 0
  )
  testthat::expect_true(
    max(bc.object.down@data, na.rm = TRUE) <= 0
  )
  ### Check that the slot @switch.point is a numeric vector.
  testthat::expect_equal(
    class(bc.object@switch.point),
    "numeric"
  )
  ### Check that the @switch.point vector doesn't contain NA values.
  testthat::expect_true(
    all(!is.na(bc.object@switch.point))
  )
  ### Check that the values in the slot @switch.point range between 0 and 1.
  testthat::expect_true(
    min(bc.object@switch.point) >= 0
  )
  testthat::expect_true(
    max(bc.object@switch.point) <= 1
  )
  ### Check that the slot @ranks is empty.
  testthat::expect_equal(
    bc.object@ranks,
    list()
  )
  ### Check that the slot @expr.matrix contains the input matrix.
  testthat::expect_equal(
    bc.object@expr.matrix,
    mtx
  )
  testthat::expect_equal(
    bc.object.up@expr.matrix,
    mtx
  )
  testthat::expect_equal(
    bc.object.down@expr.matrix,
    mtx
  )
  # Check that the slot @meta.data contains the input metadata.
  metadata <- pbmc@meta.data
  testthat::expect_equal(
    bc.object@meta.data,
    metadata
  )
  testthat::expect_equal(
    bc.object.up@meta.data,
    metadata
  )
  testthat::expect_equal(
    bc.object.down@meta.data,
    metadata
  )
  ### @SeuratInfo
  ### Check that the slot @background is empty.
  testthat::expect_equal(
    bc.object@background,
    as.matrix(data.frame())
  )
  ### Check that the slot @reductions is empty.
  testthat::expect_equal(
    bc.object@reductions,
    list()
  )
  ### Check that the slot @regression is empty.
  testthat::expect_equal(
    bc.object@regression,
    list(order = rep("", 2), vars = NULL, order.background = rep("", 2))
  )
  ### Check that the slot @n.genes is equal to the slot @n.genes in the input 
  ### geneset.
  testthat::expect_equal(
    bc.object@n.genes,
    gs10@n.genes
  )
  testthat::expect_equal(
    bc.object.up@n.genes,
    gs10up@n.genes
  )
  testthat::expect_equal(
    bc.object.down@n.genes,
    gs10down@n.genes
  )
  ### Check that the slot @mode is equal to the slot @mode in the input geneset.
  testthat::expect_equal(
    bc.object@mode,
    gs10@mode
  )
  testthat::expect_equal(
    bc.object.up@mode,
    gs10up@mode
  )
  testthat::expect_equal(
    bc.object.down@mode,
    gs10down@mode
  )
  ### Check the value of the slot @thres.
  testthat::expect_equal(
    unique(bc.object@thres, bc.object.up@thres, bc.object.down@thres),
    0.1
  )
})