test_that("default registry has the PPIs of the study and is valid", {
  reg <- default_ppi_registry()
  expect_setequal(names(reg), c("CV_t", "RN", "RNN", "D_slope", "RC", "gPi", "PPF", "PPi", "PImd", "PILSM",
                                "RTR", "PIR", "RDPI", "ESPI", "ESPIID", "PSI", "RPI", "PQ", "PR", "NRW",
                                "ESP", "CEV", "PFI", "APC", "SI", "RSI", "EVS", "FW"))
  expect_true(validate_ppi_registry(reg))
  expect_equal(reg$FW$level, "population")
})

test_that("a new PPI is added with one call", {
  my_range <- function(trait_values) max(trait_values) - min(trait_values)
  reg <- register_ppi(default_ppi_registry(), ppi_entry("MyRange", my_range))
  expect_true("MyRange" %in% names(reg))
})

test_that("wrong argument names or inputs are caught when registering", {
  f <- function(trait_values, env_values) 1
  expect_error(register_ppi(list(), ppi_entry("X", f, inputs = c(trait_values = "trait", env = "env"))),
               "function has no argument env")
  expect_error(register_ppi(list(), ppi_entry("X", f, inputs = c(trait_values = "phenotype"))), "unknown input phenotype")
})
