# Maximum n.genes in a collection.
n.max <- 500

# Filters names.
filters.names <- c("drugs", "IDs", "MoAs", "targets", "studies")

# Tests.
testthat::test_that("errors", {
  ### Check x.
  testthat::expect_error(
    GetCollection("a"),
    'x must be either PSc, SSc or DSS.'
  )
  ### Check n.genes.
  testthat::expect_error(
    GetCollection(SSc, n.genes = "a"),
    'n.genes must be numeric.'
  )
  testthat::expect_error(
    GetCollection(SSc, n.genes = 1:2),
    'n.genes must be a positive integer.'
  )
  testthat::expect_error(
    GetCollection(SSc, n.genes = 1.5),
    'n.genes must be a positive integer.'
  )
  testthat::expect_error(
    GetCollection(SSc, n.genes = -2),
    'n.genes must be a positive integer.'
  )
  testthat::expect_error(
    GetCollection(SSc, n.genes = 1000),
    paste0('n.genes exceeds the maximum number of genes in signature (',
           n.max, ').'), fixed = TRUE
  )
  ### Check mode.
  testthat::expect_error(
    GetCollection(SSc, mode = "mid"),
    'Incorrect mode.'
  )
  testthat::expect_error(
    GetCollection(SSc, mode = 1),
    'Incorrect mode.'
  )
  ### Check filters.
  testthat::expect_error(
    GetCollection(SSc, filters = "ERBB2"),
    'filters must be a list.'
  )
  testthat::expect_error(
    GetCollection(SSc, filters = list(targets = "ERBB2", foo = "foo")),
    'Invalid names in filters.'
  )
  testthat::expect_error(
    GetCollection(SSc, filters = list(targets = "ERBB2", 
                                      IDs = data.frame(0), MoAs = 1)),
    paste('Incorrect value for filter\'s entry: IDs, MoAs. You must', 
          'provide a character vector.')
  )
  ### Check include.pathways.
  testthat::expect_error(
    GetCollection(SSc, include.pathways = rep(TRUE, 2)),
    'include.pathways must be TRUE or FALSE.'
  )
  testthat::expect_error(
    GetCollection(SSc, include.pathways = 1),
    'include.pathways must be TRUE or FALSE.'
  )
})
