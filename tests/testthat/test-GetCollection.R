# --- Functions ---
# Function to aggregate and order data.frame in the slot @info.
aggregateInfo <- function(x) {
  x <- aggregate(.~ IDs, data = x, na.action = NULL, FUN = function(rw) {
    paste(na.omit(unique(rw)), collapse = ", ")
  })
  x <- x[order(x$IDs, decreasing = FALSE), ]
  return(x)
}

# --- Code ---
# Maximum n.genes in a collection.
n.max <- 250

# Filters names.
filters.names <- c("drugs", "IDs", "MoAs", "targets", "studies")

# Pathways
load("../../R/sysdata.rda")

# Test errors.
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
    GetCollection(SSc, n.genes = 500),
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

# Test warnings.
testthat::test_that("warnings", {
  ### Check when there aren't filters that yield to any result.
  testthat::expect_warning(
    GetCollection(SSc, filters = list(targets = "ERBB2", studies = "CCLE")),
    paste('The following filters\' values yielded no results:\n', 
          '- studies: CCLE\n')
  )
  testthat::expect_warning(
    GetCollection(SSc, filters = list(targets = "ERBB2", 
                                      studies = c("CCLE", "PRISM"))),
    paste('The following filters\' values yielded no results:\n', 
          '- studies: CCLE\n')
  )
  testthat::expect_warning(
    GetCollection(SSc, filters = list(targets = "DNApol", studies = "CCLE", 
                                      IDs = "sig-20879")),
    paste('The following filters\' values yielded no results:\n', 
          '- targets: DNApol\n',
          '- studies: CCLE\n')
  )
})

# Test values.
testthat::test_that("default values", {
  ### Check that GetCollection's output is a geneset object.
  testthat::expect_s4_class(
    GetCollection(SSc),
    "geneset"
  )
  ### Check that default GetCollection's output is identical to pre-loaded 
  ### collections.
  testthat::expect_equal(
    GetCollection(PSc, n.genes = 250, include.pathways = FALSE),
    PSc
  )
  testthat::expect_equal(
    GetCollection(SSc, n.genes = 250, include.pathways = FALSE),
    SSc
  )
  testthat::expect_equal(
    GetCollection(DSS, n.genes = 250, include.pathways = FALSE),
    DSS
  )
  ### Check that default GetCollection's pathways output is identical to 
  ### pathways.
  testthat::expect_equal(
    GetCollection(SSc, include.pathways = TRUE)@genelist[names(pathways)],
    pathways
  )
  ### Check the value of the slot @n.genes.
  testthat::expect_equal(
    GetCollection(SSc, n.genes = 100)@n.genes,
    100
  )
  ### Check that each character vector inside the slot @genelist has a max
  ### length of n.genes (only applies to drug signatures).
  testthat::expect_true(
    all(sapply(GetCollection(SSc, n.genes = 100, 
                             include.pathways = FALSE)@genelist, 
               FUN = function(x) all(sapply(x, length) <= 100)))
  )
  ### Check that pathway signatures are not affected by the argument n.genes.
  testthat::expect_equal(
    GetCollection(SSc, n.genes = 100, 
                  include.pathways = TRUE)@genelist[names(pathways)],
    pathways
  )
  ### Check the value of the slot @mode.
  testthat::expect_equal(
    GetCollection(SSc)@mode,
    c("up", "down")
  )
  testthat::expect_equal(
    GetCollection(SSc, mode = "up")@mode,
    "up"
  )
  testthat::expect_equal(
    GetCollection(SSc, mode = "down")@mode,
    "down"
  )
  ### Check that the names of the character vectors inside the slot @genelist 
  ### are equal to mode (applies to pathway signatures).
  testthat::expect_true(
    all(sapply(GetCollection(SSc, mode = "up", 
                             include.pathways = TRUE)@genelist,
               FUN = function(x) all(names(x) == "up")))
  )
  testthat::expect_true(
    all(sapply(GetCollection(SSc, mode = "down", 
                             include.pathways = TRUE)@genelist,
               FUN = function(x) all(names(x) == "down")))
  )
  ### If the mode = c("up", "down"), pathway signatures can have only 1 mode.
  testthat::expect_true(
    all(sapply(GetCollection(SSc, mode = c("up", "down"), 
                             include.pathways = FALSE)@genelist,
               FUN = function(x) all(identical(names(x), c("up", "down")))))
  )
  testthat::expect_true(
    all(sapply(GetCollection(SSc, mode = c("up", "down"), 
                             include.pathways = TRUE)@genelist[names(pathways)],
               FUN = function(x) all(names(x) %in% c("up", "down"))))
  )
  ### Check drugs filter.
  bortezomib <- subset(SSc@info, subset = preferred.drug.names == "BORTEZOMIB" |
                         drugs == "BORTEZOMIB")
  bortezomib <- aggregateInfo(bortezomib)
  testthat::expect_equal(
    GetCollection(SSc, filters = list(drugs = "bortezomib"))@info,
    bortezomib
  )
  ### Check IDs filter.
  sig_21378 <- subset(SSc@info, subset = IDs == "sig-21378")
  sig_21378 <- aggregateInfo(sig_21378)
  testthat::expect_equal(
    GetCollection(SSc, filters = list(IDs = "sig-21378"))@info,
    sig_21378
  )
  ### Check MoAs filter.
  NFkB_drug_ids <- subset(SSc@info, subset = MoAs == "NFkB signaling inhibitor")
  unique_ids_to_test <- unique(NFkB_drug_ids$IDs)
  NFkB <- SSc@info[SSc@info$IDs %in% unique_ids_to_test, ]
  NFkB <- aggregateInfo(NFkB)

  testthat::expect_equal(
    GetCollection(SSc, filters = list(MoAs = "NFkB signaling inhibitor"))@info,
    NFkB
  )
  ### Check targets filter.
  has_psmb9 <- grepl(pattern = "PSMB9", x = SSc@info$targets)
  unique_ids_to_test <- unique(SSc@info[has_psmb9, ]$IDs) ## tibble stuff
  PSMB9 <- SSc@info[SSc@info$IDs %in% unique_ids_to_test, ]
  PSMB9 <- aggregateInfo(PSMB9)

  testthat::expect_equal(
    GetCollection(SSc, filters = list(targets = "PSMB9"))@info,
    PSMB9
  )
  ### Check studies filter.
  PRISM <- subset(SSc@info, subset = studies == "PRISM")
  PRISM <- aggregateInfo(PRISM)
  testthat::expect_equal(
    GetCollection(SSc, filters = list(studies = "PRISM"))@info,
    PRISM
  )
  ### Check multiple filters.
  multiple <- merge(bortezomib,
                    merge(sig_21378, 
                          merge(NFkB, 
                                merge(PSMB9, 
                                      PRISM, 
                                      all = TRUE
                                      ),
                                all = TRUE
                                ),
                          all = TRUE
                          ),
                    all = TRUE
                    )
  
  testthat::expect_equal(
    GetCollection(SSc, 
                  filters = list(drugs = "bortezomib", IDs = "sig-21378",
                                 MoAs = "NFkB signaling inhibitor",
                                 targets = "PSMB9", studies = "PRISM"))@info,
    multiple
  )
})