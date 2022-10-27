# Tests.
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
    GenerateGenesets("../testdata/duplicated_gene_sets.gmt"),
    'The GMT file contains duplicated gene set\'s: GENESET1_up.'
  )
  testthat::expect_error(
    GenerateGenesets("../testdata/incorrect_mode.gmt"),
    'All gene sets\' names in the GMT file must end in "_UP" or "_DOWN/_DN".'
  )
  ### Check perform.reversal.
  testthat::expect_error(
    GenerateGenesets("../testdata/correct.gmt", 
                     perform.reversal = rep(TRUE, 2)),
    'perform.reversal must be TRUE or FALSE.'
  )
  testthat::expect_error(
    GenerateGenesets("../testdata/correct.gmt", perform.reversal = 1),
    'perform.reversal must be TRUE or FALSE.'
  )
})