# Test errors.
testthat::test_that("errors", {
  ### Check x.
  testthat::expect_error(
    GenerateGenesets("a"),
    'Pathway information must be a .gmt file'
  )
  testthat::expect_error(
    GenerateGenesets("a.gmt"),
    'cannot open file \'a.gmt\': No such file or directory'
  )
  testthat::expect_error(
    GenerateGenesets("../testdata/duplicated10.gmt"),
    paste('The GMT file contains duplicated gene set\'s: sig-20965_up,', 
          'sig-20965_down.')
  )
  testthat::expect_error(
    GenerateGenesets("../testdata/incorrect_mode10.gmt"),
    'All gene sets\' names in the GMT file must end in "_UP" or "_DOWN/_DN".'
  )
  ### Check perform.reversal.
  testthat::expect_error(
    GenerateGenesets("../testdata/correct10.gmt", 
                     perform.reversal = rep(TRUE, 2)),
    'perform.reversal must be TRUE or FALSE.'
  )
  testthat::expect_error(
    GenerateGenesets("../testdata/correct10.gmt", perform.reversal = 1),
    'perform.reversal must be TRUE or FALSE.'
  )
})

# Test values.
testthat::test_that("default values", {
  ### Test that GenerateGenesets' output is a geneset object.
  testthat::expect_s4_class(
    GenerateGenesets("../testdata/correct10.gmt"),
    "geneset"
  )
  ### Check that the slot @genelist has names.
  gs_genelist <- GenerateGenesets("../testdata/correct10.gmt")@genelist
  testthat::expect_false(
    is.null(names(gs_genelist)),
    NULL
  )
  ### Check that the slot @genelist is a list of lists.
  testthat::expect_true(
    all(sapply(gs_genelist, class) == "list")
  )
  ### Check that each list inside the slot @genelist has length 1 or 2.
  testthat::expect_true(
    all(sapply(gs_genelist, length) %in% 1:2)
  )
  ### Check that the slot @genelist is a list of lists of character vectors.
  testthat::expect_true(
    all(sapply(gs_genelist, 
               FUN = function(x) all(sapply(x, class) == "character")))
  )
  ### Check that each character vector inside the slot @genelist has a min 
  ### length of 1.
  testthat::expect_true(
    all(sapply(gs_genelist, FUN = function(x) all(sapply(x, length) >= 1)))
  )
  ### Check that each character vector inside the slot @genelist is named either
  ### "up" or "down".
  testthat::expect_true(
    all(sapply(gs_genelist, 
               FUN = function(x) all(names(x) %in% c("up", "down"))))
  )
  ### Check that the slot @n.genes is empty.
  testthat::expect_equal(
    GenerateGenesets("../testdata/correct10.gmt")@n.genes,
    NaN
  )
  ### Check that the slot @mode is c("up", "down").
  testthat::expect_equal(
    GenerateGenesets("../testdata/correct10.gmt")@mode,
    c("up", "down")
  )
  ### Check that the slot @info is empty.
  testthat::expect_equal(
    GenerateGenesets("../testdata/correct10.gmt")@info,
    data.frame()
  )
  ### Check the slot @inverse.score.
  testthat::expect_false(
    GenerateGenesets("../testdata/correct10.gmt")@inverse.score
  )
  testthat::expect_true(
    GenerateGenesets("../testdata/correct10.gmt", 
                     perform.reversal = TRUE)@inverse.score
  )
})