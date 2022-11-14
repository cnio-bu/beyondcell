# GMT filepaths.
correct.gmt.path <- "../testdata/gmt/correct10.gmt"
correct.gmt.up.path <- "../testdata/gmt/correct10up.gmt"
correct.gmt.down.path <- "../testdata/gmt/correct10down.gmt"
duplicated.gmt.path <- "../testdata/gmt/duplicated10.gmt"
incorrect.mode.gmt.path <- "../testdata/gmt/incorrect_mode10.gmt"

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
    GenerateGenesets(duplicated.gmt.path),
    paste('The GMT file contains duplicated gene set\'s: sig-20988_up,', 
          'sig-20988_down.')
  )
  testthat::expect_error(
    GenerateGenesets(incorrect.mode.gmt.path),
    'All gene sets\' names in the GMT file must end in "_UP" or "_DOWN/_DN".'
  )
  ### Check perform.reversal.
  testthat::expect_error(
    GenerateGenesets(correct.gmt.path, perform.reversal = rep(TRUE, 2)),
    'perform.reversal must be TRUE or FALSE.'
  )
  testthat::expect_error(
    GenerateGenesets(correct.gmt.path, perform.reversal = 1),
    'perform.reversal must be TRUE or FALSE.'
  )
})

# Test messages.
testthat::test_that("messages", {
  ### Check initial message.
  testthat::expect_message(
    GenerateGenesets(correct.gmt.path),
    'Reading gmt file...'
  )
})

# Test values.
testthat::test_that("default values", {
  ### Test that GenerateGenesets' output is a geneset object.
  testthat::expect_s4_class(
    GenerateGenesets(correct.gmt.path),
    "geneset"
  )
  ### Check that the slot @genelist has names.
  gs_genelist <- GenerateGenesets(correct.gmt.path)@genelist
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
  ### Check the values of the slot @n.genes.
  testthat::expect_equal(
    GenerateGenesets(correct.gmt.path)@n.genes,
    100 + 100
  )
  ### Check the values of the slot @mode.
  testthat::expect_equal(
    GenerateGenesets(correct.gmt.path)@mode,
    c("up", "down")
  )
  testthat::expect_equal(
    GenerateGenesets(correct.gmt.up.path)@mode,
    "up"
  )
  testthat::expect_equal(
    GenerateGenesets(correct.gmt.down.path)@mode,
    "down"
  )
  ### Check that the slot @info is empty.
  testthat::expect_equal(
    GenerateGenesets(correct.gmt.path)@info,
    data.frame()
  )
  ### Check the slot @inverse.score.
  testthat::expect_false(
    GenerateGenesets(correct.gmt.path)@inverse.score
  )
  testthat::expect_true(
    GenerateGenesets(correct.gmt.path, perform.reversal = TRUE)@inverse.score
  )
})