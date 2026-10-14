# Load all functions of the simulation, PPI and pipeline modules
for (dir in c("sim", "ppi", "pipeline")) {
  files <- sort(list.files(here::here("R", dir), pattern = "\\.R$", full.names = TRUE))
  for (f in setdiff(files, here::here("R", "pipeline", "load_pipeline.R"))) source(f)
}
