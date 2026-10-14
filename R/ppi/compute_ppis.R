# =============================================================================
# compute_ppis.R — compute the PPIs of a registry for every shape
# =============================================================================

# Genotypes x environments matrix of replicate means (rows named by genotype)
replicate_means <- function(phenotypes, environments) {
  M <- tapply(phenotypes$Trait, list(phenotypes$Genotype, phenotypes$Environment), mean)
  M <- M[unique(phenotypes$Genotype), as.character(environments), drop = FALSE]
  matrix(M, nrow(M), ncol(M), dimnames = list(rownames(M), NULL))
}

first_value <- function(x) as.numeric(if (is.list(x)) x[[1]] else x)[1]

# Value of one PPI for every genotype of one shape (named numeric vector)
compute_ppi <- function(entry, M, environments) {
  call_fun <- function(context) {
    args <- c(lapply(entry$inputs, function(key) context[[key]]), entry$args)
    # Warnings of the index functions are silenced as in the rest of the project;
    # an error gives NA for this genotype instead of stopping the run
    # (some index functions also print progress bars, which are captured)
    res <- NA_real_
    utils::capture.output(res <- tryCatch(suppressWarnings(do.call(entry$fun, args)), error = function(e) NA_real_))
    res
  }
  if (entry$level == "population") {
    res <- call_fun(list(matrix = M, env = environments))
    values <- if (is.null(entry$extract)) first_value(res) else entry$extract(res)
    return(values[rownames(M)])
  }
  extract <- if (is.null(entry$extract)) first_value else entry$extract
  context <- list(env = environments, env_index = seq_along(environments), covariate = rep(1, length(environments)))
  values <- vapply(rownames(M), function(g) {
    context$trait <- M[g, ]
    extract(call_fun(context))
  }, numeric(1))
  n_failed <- sum(is.na(values))
  if (n_failed) message("PPI ", entry$name, ": no value for ", n_failed, " genotype(s)")
  values
}

# PPI values of every genotype (rows) and replicate (columns) of one shape:
# the PPI is computed on each replicate's own reaction norm
ppi_per_replicate <- function(entry, phenotypes, environments) {
  ids <- unique(phenotypes$Genotype)
  reps <- sort(unique(phenotypes$Replicate))
  out <- sapply(reps, function(r) {
    M <- replicate_means(phenotypes[phenotypes$Replicate == r, ], environments)
    compute_ppi(entry, M, environments)[ids]
  })
  matrix(out, length(ids), length(reps), dimnames = list(ids, reps))
}

# phenotypes: named list (one table per shape). Returns a named list (one table per
# PPI) with columns Score, Type (shape), Genotype_ID (line ID) and Replicate.
# replicate_handling decides on which level a PPI value is computed:
#   "mean_before"  replicates are averaged per genotype and environment first,
#                  then one PPI value per genotype (Replicate = NA); the standard
#   "replicates"   one PPI value per replicate's reaction norm, kept as separate rows
compute_ppis <- function(phenotypes, environments, registry = default_ppi_registry(), exclude = NULL,
                         replicate_handling = c("mean_before", "replicates")) {
  replicate_handling <- match.arg(replicate_handling)
  unknown <- setdiff(exclude, names(registry))
  if (length(unknown)) stop("cannot exclude unknown PPI(s): ", paste(unknown, collapse = ", "))
  registry <- registry[setdiff(names(registry), exclude)]
  lapply(registry, function(entry) {
    do.call(rbind, lapply(names(phenotypes), function(shape) {
      ph <- phenotypes[[shape]]
      if (replicate_handling == "mean_before") {
        M <- replicate_means(ph, environments)
        v <- compute_ppi(entry, M, environments)
        return(data.frame(Score = unname(v), Type = shape, Genotype_ID = rownames(M), Replicate = NA_integer_,
                          stringsAsFactors = FALSE))
      }
      P <- ppi_per_replicate(entry, ph, environments)
      data.frame(Score = as.vector(t(P)), Type = shape, Genotype_ID = rep(rownames(P), each = ncol(P)),
                 Replicate = rep(as.integer(colnames(P)), times = nrow(P)), stringsAsFactors = FALSE)
    }))
  })
}
