# Run all unit tests. Usage (from the repo root): Rscript --vanilla tests/run_tests.R
library(testthat)
test_dir("tests/testthat", reporter = "summary", stop_on_failure = TRUE)
