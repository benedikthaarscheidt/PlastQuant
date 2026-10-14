# Load every module file once per test run
root <- normalizePath(file.path(testthat::test_path(), "..", ".."))
for (dir in c("sim", "ppi", "pipeline")) {
  files <- sort(list.files(file.path(root, "R", dir), pattern = "\\.R$", full.names = TRUE))
  for (f in files[basename(files) != "load_pipeline.R"]) source(f, local = FALSE)
}
