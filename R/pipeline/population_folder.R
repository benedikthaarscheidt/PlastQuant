# =============================================================================
# population_folder.R — one folder per population: its data and everything
#                       computed from it
# =============================================================================
# <population>/
#   data/            the population (simulated or own data)
#     traits/<trait>.csv  first column Genotype (one row per plant: a genotype appears
#                         once per replicate), then one column per environment
#     environments.csv    Column, Value: numeric environment value of every trait column
#     snp_matrix.csv      Genotype, then one column per SNP, dosages 0/1/2 (optional;
#                         simulated populations store snp_matrix.rds)
#     kinship.csv         Genotype, then one column per genotype (optional; otherwise
#                         computed from the SNPs)
#     truth.csv           Trait, Parameter, SNP: known causal SNPs (optional)
#     config.rds, mother_parameters.csv, simulation.rds   simulated populations only
#   ppis/  gwas/  variance/  analyses/<analysis>/  plots/   results (rewritten by each run)
# =============================================================================

# Subfolder of a population folder, created if needed
population_subdir <- function(dir, ...) {
  d <- file.path(dir, ...)
  dir.create(d, recursive = TRUE, showWarnings = FALSE)
  d
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

# Write a population (output of simulate_reaction_norms() or pipeline_input()) into
# dir/data/. Returns dir.
write_population <- function(x, dir, snp_format = if (is.null(x$config)) "csv" else "rds") {
  data_dir <- population_subdir(dir, "data")
  unlink(file.path(data_dir, "traits"), recursive = TRUE)
  traits_dir <- population_subdir(dir, "data", "traits")
  env_map <- setNames(x$environments, paste0("E", seq_along(x$environments)))
  write.csv(data.frame(Column = names(env_map), Value = unname(env_map)),
            file.path(data_dir, "environments.csv"), row.names = FALSE)
  for (trait in names(x$phenotypes)) {
    write.csv(long_to_trait_matrix(x$phenotypes[[trait]], env_map),
              file.path(traits_dir, paste0(trait, ".csv")), row.names = FALSE)
  }
  if (!is.null(x$snp_matrix)) {
    if (snp_format == "rds") {
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
  if (!is.null(x$config)) {
    saveRDS(x$config, file.path(data_dir, "config.rds"))
    write.csv(x$mother_parameters, file.path(data_dir, "mother_parameters.csv"), row.names = FALSE)
    saveRDS(x, file.path(data_dir, "simulation.rds"))
  }
  dir
}

read_genotype_table <- function(path) {
  tab <- read.csv(path, check.names = FALSE, stringsAsFactors = FALSE, colClasses = c(Genotype = "character"))
  if (!"Genotype" %in% names(tab)) stop(path, ": the first column must be Genotype")
  tab
}

# Read dir/data/ into a pipeline input (see pipeline_input.R) plus the kinship matrix
# and, for simulated populations, config and mother_parameters.
read_population <- function(dir) {
  data_dir <- file.path(dir, "data")
  trait_files <- sort(list.files(file.path(data_dir, "traits"), pattern = "\\.csv$", full.names = TRUE))
  if (!length(trait_files)) stop("no data/traits/<trait>.csv files in ", dir)
  env_file <- file.path(data_dir, "environments.csv")
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
  if (file.exists(file.path(data_dir, "snp_matrix.rds"))) {
    snp_matrix <- readRDS(file.path(data_dir, "snp_matrix.rds"))
  } else if (file.exists(file.path(data_dir, "snp_matrix.csv"))) {
    tab <- read_genotype_table(file.path(data_dir, "snp_matrix.csv"))
    snp_matrix <- as.matrix(tab[, -1, drop = FALSE])
    rownames(snp_matrix) <- tab$Genotype
  }
  kinship <- NULL
  if (file.exists(file.path(data_dir, "kinship.csv"))) {
    tab <- read_genotype_table(file.path(data_dir, "kinship.csv"))
    kinship <- as.matrix(tab[, -1, drop = FALSE])
    rownames(kinship) <- tab$Genotype
    if (!identical(sort(rownames(kinship)), sort(colnames(kinship)))) stop("kinship.csv: row and column genotypes differ")
  }
  truth <- NULL
  if (file.exists(file.path(data_dir, "truth.csv"))) {
    truth <- read.csv(file.path(data_dir, "truth.csv"), stringsAsFactors = FALSE)
    names(truth)[names(truth) == "Trait"] <- "Shape"
  }

  pop <- pipeline_input(phenotypes, sort(unique(unname(env_map))), snp_matrix, truth)
  pop$kinship <- kinship
  if (!is.null(kinship)) {
    missing <- setdiff(unique(unlist(lapply(phenotypes, `[[`, "Genotype"))), rownames(kinship))
    if (length(missing)) stop("kinship.csv lacks genotypes: ", paste(head(missing, 5), collapse = ", "))
  }
  if (file.exists(file.path(data_dir, "config.rds"))) {
    pop$config <- readRDS(file.path(data_dir, "config.rds"))
    pop$mother_parameters <- read.csv(file.path(data_dir, "mother_parameters.csv"), stringsAsFactors = FALSE)
    # simulated shapes in the order of the model, not alphabetically
    pop$phenotypes <- pop$phenotypes[intersect(pop$config$shapes, names(pop$phenotypes))]
  }
  pop$name <- basename(normalizePath(dir))
  pop
}

# Stop unless the population in dir has at least n measurement points (for analyses
# that subsample a dense curve)
require_population_points <- function(dir, n) {
  env_file <- file.path(dir, "data", "environments.csv")
  if (!file.exists(env_file)) stop("no population in ", dir, " (simulate_population() or place own data in data/)")
  have <- nrow(read.csv(env_file))
  if (have < n) stop("this analysis needs a population with at least ", n, " measurement points; ",
                     dir, " has ", have, " (simulate one with sim_config(..., n_points = ", n, "))")
  invisible(have)
}

# The reaction norms of a simulated population, sampled again at n_points equidistant
# points (same genotypes and plants; see measure_norms()). Own data cannot be
# resampled: its reaction norms are only known at the measured environments.
sample_population <- function(dir, n_points, error = TRUE) {
  sim_file <- file.path(dir, "data", "simulation.rds")
  if (!file.exists(sim_file)) stop(dir, " is not a simulated population (no data/simulation.rds); own data cannot be resampled")
  measure_norms(readRDS(sim_file), n_points, error)
}
