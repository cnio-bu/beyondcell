expected_gene_set_names <- scan("../testdata/gmt/correct10_names.txt",
                                what = character())

total_genesets <- length(expected_gene_set_names)
expected_genes_geneset <- rep(100, times = total_genesets)
names(expected_genes_geneset) <- expected_gene_set_names

# Test that readGMT returns a list of gene sets
test_that("readGMT returns a list", {
  result <- readGMT("../testdata/gmt/correct10.gmt")
  expect_true(is.list(result))
})

# Test that readGMT returns the correct number of gene sets
# Remember that genesets are bidirectional so correct10 = 20 genesets (UP/DN)
test_that("readGMT returns the correct number of gene sets", {
  result <- readGMT("../testdata/gmt/correct10.gmt")
  expect_equal(length(result), total_genesets)
})

# Test that readGMT returns gene sets with the correct names
test_that("readGMT returns gene sets with the correct names", {
  result <- readGMT("../testdata/gmt/correct10.gmt")
  expect_equal(names(result), expected_gene_set_names)
})

# Test that readGMT returns gene sets with the correct number of genes
test_that("readGMT returns gene sets with the correct number of genes", {
  result <- readGMT("../testdata/gmt/correct10.gmt")
  result_lengths <- sapply(result, length)
  expect_equal(result_lengths, expected_genes_geneset)
})

# Test that readGMT handles invalid input correctly
test_that("readGMT handles invalid input correctly", {
  expect_error(readGMT(123), "x must be a single string.")
  expect_error(readGMT("../testdata/gmt/unknown.gmt"), "../testdata/gmt/unknown.gmt does not exist.")
  expect_error(readGMT("../testdata/single-cell/genes.tsv"), "../testdata/single-cell/genes.tsv must be a GMT file.")
})
