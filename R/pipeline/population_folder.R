# =============================================================================
# population_folder.R — one folder per population: the population, its sampling
#                       at different resolutions, and everything computed from them
# =============================================================================
# <population>/
#   data/              the population itself (no measurements)
#     norms.rds              true mother and replicate reaction norms (simulated only)
#     mother_parameters.csv, plant_parameters.csv, config.rds          (simulated only)
#     snp_matrix.csv         Genotype, then one column per SNP, dosages 0/1/2 (optional;
#                            simulated populations store snp_matrix.rds)
#     kinship.csv            Genotype, then one column per genotype (optional; otherwise
#                            computed from the SNPs)
#     truth.csv              Trait, Parameter, SNP: known causal SNPs (optional)
#   variance/          latent variance components (simulated only)
#   n_<k>/             one resolution: the population measured at k environments
#     data/traits/<trait>.csv  first column Genotype (one row per plant: a genotype
#                              appears once per replicate), then one column per environment
#     data/environments.csv    Column, Value: numeric environment value of every column
#     variance/  ppis/  gwas/  analyses/  plots/      results at this resolution
#   analyses/  plots/  results across resolutions
# =============================================================================

# Subfolder, created if needed
population_subdir <- function(dir, ...) {
  d <- file.path(dir, ...)
  dir.create(d, recursive = TRUE, showWarnings = FALSE)
  d
}

# Folder of the resolution with n_points sampling points
resolution_path <- function(population_dir, n_points) file.path(population_dir, paste0("n_", n_points))

# The resolution folders of a population (n_<k> with measurements), ordered by k
resolution_dirs <- function(population_dir) {
  d <- list.dirs(population_dir, full.names = TRUE, recursive = FALSE)
  d <- d[grepl("^n_[0-9]+$", basename(d)) & dir.exists(file.path(d, "data", "traits"))]
  d[order(as.integer(sub("^n_", "", basename(d))))]
}

# Trait matrix (data.frame: Genotype, then one column per environment) -> long table.
# env_map: named numeric vector, names = trait columns, values = environment values.
# Replicates are numbered in row order within each genotype.
trait_matrix_to_long <- function(tab, env_map) {
  ids <- as.character(tab$Genotype)
  rep_no <- as.integer(ave(seq_along(ids), ids, FUN = seq_along))
  cols <- names(env_map)
  data.frame(
    Genotype    = rep(ids, each = length(cols)),
    Replicate   = rep(rep_no, each = length(cols)),
    Environment = rep(unname(env_map), times = length(ids)),
    Trait       = as.vector(t(as.matrix(tab[, cols, drop = FALSE]))),
    stringsAsFactors = FALSE
  )
}

# Long table -> trait matrix (one row per genotype and replicate, columns names(env_map))
long_to_trait_matrix <- function(ph, env_map) {
  ph <- ph[order(match(ph$Genotype, unique(ph$Genotype)), ph$Replicate), ]
  plants <- unique(ph[, c("Genotype", "Replicate")])
  key <- paste(ph$Genotype, ph$Replicate, sep = "\r")
  out <- data.frame(Genotype = as.character(plants$Genotype), stringsAsFactors = FALSE)
  for (col in names(env_map)) {
    sel <- ph$Environment == env_map[[col]]
    out[[col]] <- ph$Trait[sel][match(paste(plants$Genotype, plants$Replicate, sep = "\r"), key[sel])]
  }
  out
}

# Write the population into dir/data/: for a simulation (simulate_norms(), or a measured
# simulation) the true norms, parameters, SNPs, truth and settings; for own data built in
# R (pipeline_input()) the SNPs, kinship and truth, plus its measurements as the
# resolution folder n_<number of environments>. Returns dir.
write_population <- function(x, dir) {
  data_dir <- population_subdir(dir, "data")
  simulated <- !is.null(x$config)
  if (!is.null(x$snp_matrix)) {
    if (simulated) {
      saveRDS(x$snp_matrix, file.path(data_dir, "snp_matrix.rds"))
    } else {
      write.csv(data.frame(Genotype = rownames(x$snp_matrix), x$snp_matrix, check.names = FALSE),
                file.path(data_dir, "snp_matrix.csv"), row.names = FALSE)
    }
  }
  if (!is.null(x$kinship)) {
    write.csv(data.frame(Genotype = rownames(x$kinship), x$kinship, check.names = FALSE),
              file.path(data_dir, "kinship.csv"), row.names = FALSE)
  }
  if (!is.null(x$truth)) {
    tr <- x$truth
    names(tr)[names(tr) == "Shape"] <- "Trait"
    write.csv(tr, file.path(data_dir, "truth.csv"), row.names = FALSE)
  }
  if (simulated) {
    norms <- x[setdiff(names(x), c("environments", "phenotypes", "mother_curves", "measurement_sd"))]
    saveRDS(norms, file.path(data_dir, "norms.rds"))
    saveRDS(x$config, file.path(data_dir, "config.rds"))
    write.csv(x$mother_parameters, file.path(data_dir, "mother_parameters.csv"), row.names = FALSE)
    # one table for all shapes; parameters a shape does not have are NA
    cols <- unique(unlist(lapply(x$plant_parameters, names)))
    plants <- do.call(rbind, lapply(names(x$plant_parameters), function(s) {
      d <- x$plant_parameters[[s]]
      for (col in setdiff(cols, names(d))) d[[col]] <- NA_real_
      cbind(Shape = s, d[, cols])
    }))
    write.csv(plants, file.path(data_dir, "plant_parameters.csv"), row.names = FALSE)
  } else if (!is.null(x$phenotypes)) {
    write_resolution(x, resolution_path(dir, length(x$environments)))
  }
  dir
}

# Write one measurement (measure_norms() / sample_population(), or own phenotypes and
# environments) into the resolution folder res_dir/data/. Returns res_dir.
write_resolution <- function(measured, res_dir) {
  unlink(file.path(res_dir, "data"), recursive = TRUE)
  traits_dir <- population_subdir(res_dir, "data", "traits")
  env_map <- setNames(measured$environments, paste0("E", seq_along(measured$environments)))
  write.csv(data.frame(Column = names(env_map), Value = unname(env_map)),
            file.path(res_dir, "data", "environments.csv"), row.names = FALSE)
  for (trait in names(measured$phenotypes)) {
    write.csv(long_to_trait_matrix(measured$phenotypes[[trait]], env_map),
              file.path(traits_dir, paste0(trait, ".csv")), row.names = FALSE)
  }
  res_dir
}

read_genotype_table <- function(path) {
  tab <- read.csv(path, check.names = FALSE, stringsAsFactors = FALSE, colClasses = c(Genotype = "character"))
  if (!"Genotype" %in% names(tab)) stop(path, ": the first column must be Genotype")
  tab
}

# Read a resolution folder (<population>/n_<k>): its measurements, plus the SNPs,
# kinship, truth and (simulated) settings and parameters of the population. Returns a
# pipeline input (see pipeline_input.R) with kinship, config, mother_parameters,
# population_dir and resolution_dir.
read_population <- function(res_dir) {
  res_data <- file.path(res_dir, "data")
  if (!grepl("^n_[0-9]+$", basename(res_dir)) || !dir.exists(res_data)) {
    stop(res_dir, " is not a resolution folder (<population>/n_<points>/data; create one with sample_resolution())")
  }
  pop_dir <- dirname(normalizePath(res_dir))
  pop_data <- file.path(pop_dir, "data")

  trait_files <- sort(list.files(file.path(res_data, "traits"), pattern = "\\.csv$", full.names = TRUE))
  if (!length(trait_files)) stop("no data/traits/<trait>.csv files in ", res_dir)
  env_file <- file.path(res_data, "environments.csv")
  if (!file.exists(env_file)) stop("missing ", env_file, " (columns Column, Value)")
  env_tab <- read.csv(env_file, stringsAsFactors = FALSE)
  env_map <- setNames(as.numeric(env_tab$Value), env_tab$Column)

  phenotypes <- lapply(trait_files, function(f) {
    first <- names(read.csv(f, nrows = 1, check.names = FALSE))[1]
    if (first != "Genotype") stop(basename(f), ": the first column must be Genotype (the genotype ID of each row)")
    tab <- read_genotype_table(f)
    unknown <- setdiff(names(tab)[-1], names(env_map))
    if (length(unknown)) stop(basename(f), ": columns not listed in environments.csv: ", paste(unknown, collapse = ", "))
    trait_matrix_to_long(tab, env_map[names(tab)[-1]])
  })
  names(phenotypes) <- sub("\\.csv$", "", basename(trait_files))

  snp_matrix <- NULL
  if (file.exists(file.path(pop_data, "snp_matrix.rds"))) {
    snp_matrix <- readRDS(file.path(pop_data, "snp_matrix.rds"))
  } else if (file.exists(file.path(pop_data, "snp_matrix.csv"))) {
    tab <- read_genotype_table(file.path(pop_data, "snp_matrix.csv"))
    snp_matrix <- as.matrix(tab[, -1, drop = FALSE])
    rownames(snp_matrix) <- tab$Genotype
  }
  kinship <- NULL
  if (file.exists(file.path(pop_data, "kinship.csv"))) {
    tab <- read_genotype_table(file.path(pop_data, "kinship.csv"))
    kinship <- as.matrix(tab[, -1, drop = FALSE])
    rownames(kinship) <- tab$Genotype
    if (!identical(sort(rownames(kinship)), sort(colnames(kinship)))) stop("kinship.csv: row and column genotypes differ")
  }
  truth <- NULL
  if (file.exists(file.path(pop_data, "truth.csv"))) {
    truth <- read.csv(file.path(pop_data, "truth.csv"), stringsAsFactors = FALSE)
    names(truth)[names(truth) == "Trait"] <- "Shape"
  }

  pop <- pipeline_input(phenotypes, sort(unique(unname(env_map))), snp_matrix, truth)
  pop$kinship <- kinship
  if (!is.null(kinship)) {
    missing <- setdiff(unique(unlist(lapply(phenotypes, `[[`, "Genotype"))), rownames(kinship))
    if (length(missing)) stop("kinship.csv lacks genotypes: ", paste(head(missing, 5), collapse = ", "))
  }
  if (file.exists(file.path(pop_data, "config.rds"))) {
    pop$config <- readRDS(file.path(pop_data, "config.rds"))
    pop$mother_parameters <- read.csv(file.path(pop_data, "mother_parameters.csv"), stringsAsFactors = FALSE)
    # simulated shapes in the order of the model, not alphabetically
    pop$phenotypes <- pop$phenotypes[intersect(pop$config$shapes, names(pop$phenotypes))]
  }
  pop$population_dir <- pop_dir
  pop$resolution_dir <- normalizePath(res_dir)
  pop
}

# The true reaction norms of a simulated population (data/norms.rds) sampled at
# n_points equidistant points or at given points t in [0, 1] (same genotypes and
# plants; see measure_norms()). Own data cannot be resampled: its reaction norms are
# only known at the measured environments.
sample_population <- function(dir, n_points = NULL, error = TRUE, t = NULL) {
  norms_file <- file.path(dir, "data", "norms.rds")
  if (!file.exists(norms_file)) stop(dir, " is not a simulated population (no data/norms.rds); own data cannot be resampled")
  measure_norms(readRDS(norms_file), n_points, error, t)
}
