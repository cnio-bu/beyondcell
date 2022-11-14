# --- Functions ---
# Function to capitalize genes (mouse symbol-like).
capitalize <- function(x) {
  x <- tolower(x)
  substr(x, 1, 1) <- toupper(substr(x, 1, 1))
  return(x)
}

# Function to compute the proportion of expressed genes per cell and signature.
proportion.expressed <- function(mtx, genelist) {
  sapply(genelist, FUN = function(x) {
    genes <- intersect(rownames(mtx), x)
    apply(mtx[genes, ], 2, function(y) sum(y > 0)/length(x))
  })
}

# Function to add NaN to proportions below an expression threshold.
below.thres.nan <- function(prop, expr.thres) {
  prop[prop < expr.thres] <- NaN
  return(prop)
}

# Function to compute the normalized BCS for a certain mode.
normalized.score <- function(mtx, gs, mode) {
  gs <- lapply(gs@genelist, FUN = function(x) unique(x[[mode]]))
  t(sapply(gs, FUN = function(y) {
    genes <- intersect(rownames(mtx), y)
    raw <- colMeans(mtx[genes, ])
    sum.expr <- colSums(mtx[genes, ])
    sig.stdev <- apply(mtx[genes, ], 2, sd)
    raw * ((sum.expr - sig.stdev)/(raw + sig.stdev))
  }))
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

# Geneset with only one signature.
gs1 <- GetCollection(SSc, n.genes = 100, mode = c("up", "down"),
                     filters = list(IDs = SSc@info$IDs[1]),
                     include.pathways = FALSE)

# Pathways
load("../../R/sysdata.rda")

# Visium data with SCT normalization.
image <- Seurat::Read10X_Image(image.dir = "../testdata/visium/", 
                               image.name = "tissue_lowres_image.png")
pbmc.visium <- Seurat::Load10X_Spatial(data.dir = "../testdata/visium/", 
                                       filename = "matrix.h5", image = image)
pbmc.visium <- SCTransform(pbmc.visium, assay = "Spatial", verbose = FALSE)

# Complete cases.
gs10.genelist <- lapply(gs10@genelist, FUN = function(x) unique(unlist(x)))
proportion <- proportion.expressed(mtx, gs10.genelist)
complete.cases.01 <- sum(complete.cases(below.thres.nan(proportion, 0.1)))
complete.cases.03 <- sum(complete.cases(below.thres.nan(proportion, 0.3)))

# Normalized scores.
nan.01 <- t(proportion) < 0.1

norm.up <- normalized.score(mtx, gs10, mode = "up")
norm.up[nan.01] <- 0
norm.down <- -1 * normalized.score(mtx, gs10, mode = "down")
norm.down[nan.01] <- 0

norm <- round(norm.up + norm.down, digits = 2)
norm.up <- round(norm.up, digits = 2)
norm.down <- round(norm.down, digits = 2)
norm.up[nan.01] <- NaN
norm.down[nan.01] <- NaN
norm[nan.01 | is.na(norm.up) & is.na(norm.down)] <- NaN

# Scaled scores.
scaled <- scale.score(norm)
scaled.up <- scale.score(norm.up)
scaled.down <- scale.score(norm.down)

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
  ### Check that bcScores does not throw any error with only one signature.
  testthat::expect_no_error(
    bcScore(pbmc, gs = gs1, expr.thres = 0.1)
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
  testthat::expect_message(
    bcScore(pbmc, gs = gs10, expr.thres = 0.1),
    paste0('There are ', complete.cases.01, '/80 cells without missing values ', 
           'in your beyondcell object.')
  )
  testthat::expect_message(
    bcScore(pbmc, gs = gs10, expr.thres = 0.3),
    paste0('There are ', complete.cases.03, '/80 cells without missing values ', 
           'in your beyondcell object.')
  )
})

# Test values.
testthat::test_that("default values", {
  bc.object <- bcScore(pbmc, gs = gs10, expr.thres = 0.1)
  bc.object.up <- bcScore(pbmc, gs = gs10up, expr.thres = 0.1)
  bc.object.down <- bcScore(pbmc, gs = gs10down, expr.thres = 0.1)
  bc.object.visium <- bcScore(pbmc.visium, gs = gs10, expr.thres = 0.1)
  ### Check that bcScore output is a beyondcell object.
  testthat::expect_s4_class(
    bc.object,
    "beyondcell"
  )
  testthat::expect_s4_class(
    bc.object.visium,
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
  ### Check that the slot @normalized is equal to the slot @data.
  testthat::expect_equal(
    bc.object@normalized,
    bc.object@data
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
  ### Check that the matrices in slots @scaled, @normalized and @data have the 
  ### same dimensions.
  testthat::expect_true(
    identical(dim(bc.object@scaled), dim(bc.object@normalized))
  )
  ### Check that the matrices in slots @scaled, @normalized and @data have the 
  ### same dimnames.
  testthat::expect_true(
    identical(dimnames(bc.object@scaled), dimnames(bc.object@normalized))
  )
  ### Check that the matrices contain NaN values (not NA).
  testthat::expect_true(
    all(is.nan(bc.object@scaled[is.na(bc.object@scaled)]))
  )
  testthat::expect_true(
    all(is.nan(bc.object@normalized[is.na(bc.object@normalized)]))
  )
  ### Check that the NaN values are at the same indexes.
  testthat::expect_true(
    identical(which(is.nan(bc.object@scaled)), 
              which(is.nan(bc.object@normalized)))
  )
  ### Check that the non NaN values in the slot @scaled range between 0 and 1.
  testthat::expect_equal(
    c(min(bc.object@scaled, na.rm = TRUE), max(bc.object@scaled, na.rm = TRUE)),
    0:1
  )
  ### Check the values of the slot @scaled.
  testthat::expect_equal(
    bc.object@scaled, scaled
  )
  testthat::expect_equal(
    bc.object.up@scaled, scaled.up
  )
  testthat::expect_equal(
    bc.object.down@scaled, scaled.down
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
  ### Check the values of the slot @normalized.
  testthat::expect_equal(
    bc.object@normalized, norm
  )
  testthat::expect_equal(
    bc.object.up@normalized, norm.up
  )
  testthat::expect_equal(
    bc.object.down@normalized, norm.down
  )
  ### Check that the pathway scores are the same independently of the value of 
  ### the slot @inverse.score in the input geneset object.
  dss <- GetCollection(DSS, n.genes = 100, include.pathways = TRUE)
  ssc <- GetCollection(SSc, n.genes = 100, include.pathways = TRUE)
  bc.object.dss <- bcScore(pbmc, gs = dss, expr.thres = 0)
  bc.object.ssc <- bcScore(pbmc, gs = ssc, expr.thres = 0)
  testthat::expect_equal(
    bc.object.dss@normalized[names(pathways), ], 
    bc.object.ssc@normalized[names(pathways), ]
  )
  testthat::expect_equal(
    bc.object.dss@scaled[names(pathways), ], 
    bc.object.ssc@scaled[names(pathways), ]
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
  ### Check the values of the slot @switch.point.
  testthat::expect_equal(
    unname(bc.object@switch.point),
    sp.up.down(norm, scaled)
  )
  testthat::expect_equal(
    unique(bc.object.up@switch.point),
    0
  )
  testthat::expect_equal(
    unique(bc.object.down@switch.point),
    1
  )
  ### Check that the slot @ranks is empty.
  testthat::expect_equal(
    bc.object@ranks,
    list()
  )
  testthat::expect_equal(
    bc.object.visium@ranks,
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
  testthat::expect_equal(
    bc.object.visium@expr.matrix,
    as.matrix(pbmc.visium@assays$SCT@data)
  )
  # Check that the slot @meta.data contains the input metadata.
  metadata <- pbmc@meta.data
  testthat::expect_equal(
    bc.object@meta.data,
    metadata
  )
  testthat::expect_equal(
    bc.object.visium@meta.data,
    pbmc.visium@meta.data
  )
  ### Check the values of the slot @SeuratInfo.
  testthat::expect_equal(
    bc.object@SeuratInfo,
    list(assays = pbmc@active.assay, reductions = pbmc@reductions, 
         images = pbmc@images)
  )
  testthat::expect_equal(
    bc.object.visium@SeuratInfo,
    list(assays = pbmc.visium@active.assay, reductions = pbmc.visium@reductions, 
         images = pbmc.visium@images)
  )
  ### Check that the slot @background is empty.
  testthat::expect_equal(
    dim(bc.object@background),
    c(0, 0)
  )
  testthat::expect_equal(
    dim(bc.object.visium@background),
    c(0, 0)
  )
  ### Check that the slot @reductions is empty.
  testthat::expect_equal(
    bc.object@reductions,
    list()
  )
  testthat::expect_equal(
    bc.object.visium@reductions,
    list()
  )
  ### Check that the slot @regression is empty.
  testthat::expect_equal(
    bc.object@regression,
    list(order = rep("", 2), vars = NULL, order.background = rep("", 2))
  )
  testthat::expect_equal(
    bc.object.visium@regression,
    list(order = rep("", 2), vars = NULL, order.background = rep("", 2))
  )
  ### Check that the slot @n.genes is equal to the slot @n.genes in the input 
  ### geneset.
  testthat::expect_equal(
    bc.object@n.genes,
    gs10@n.genes
  )
  testthat::expect_equal(
    bc.object.visium@n.genes,
    gs10@n.genes
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
  ### Check the values of the slot @thres.
  testthat::expect_equal(
    unique(c(bc.object@thres, bc.object.up@thres, bc.object.down@thres, 
             bc.object.visium@thres)),
    0.1
  )
})