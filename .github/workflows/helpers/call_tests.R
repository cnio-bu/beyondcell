library("beyondcell")
library("devtools")
library("withr")

withr::with_output_sink("test_summary.txt", append = TRUE, code = {
  devtools::test(show_report = TRUE, ... = default_compact_reporter())
} )
