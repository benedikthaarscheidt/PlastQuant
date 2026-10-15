# =============================================================================
# compute_ppis.R — compute the PPIs of a registry for every shape
# =============================================================================

# Genotypes x environments matrix of replicate means (rows named by genotype)
replicate_means <- function(phenotypes, environments) {
  ids <- as.character(phenotypes$Genotype)
  M <- tapply(phenotypes$Trait, list(ids, phenotypes$Environment), mean)
  M <- M[unique(ids), as.character(environments), drop = FALSE]
  matrix(M, nrow(M), ncol(M), dimnames = list(rownames(M), NULL))
}

first_value <- function(x) as.numeric(if (is.list(x)) x[[1]] else x)[1]

# Value of one PPI for every genotype of one shape (named numeric vector). The printed
# output of the PPI functions (e.g. progress bars) is shown; warnings are silenced as in
# the rest of the project; an error is reported with the PPI, the genotype and the
# error message, and gives NA for that genotype.
compute_ppi <- function(entry, M, environments) {
  call_fun <- function(context, who) {
    args <- c(lapply(entry$inputs, function(key) context[[key]]), entry$args)
    tryCatch({
      v <- extract(suppressWarnings(do.call(entry$fun, args)))   # warnings are silenced project-wide
      setNames(as.numeric(v), names(v))
    }, error = function(e) {
      message("PPI ", entry$name, ", ", who, ": error: ", conditionMessage(e))
      NA_real_
    })
  }
  extract <- if (is.null(entry$extract)) first_value else entry$extract
  if (entry$level == "population") {
    res <- call_fun(list(matrix = M, env = environments), "all genotypes")
    values <- if (is.null(names(res))) setNames(rep(NA_real_, nrow(M)), rownames(M)) else res[rownames(M)]
    names(values) <- rownames(M)
  } else {
    context <- list(env = environments, env_index = seq_along(environments), covariate = rep(1, length(environments)))
    values <- vapply(rownames(M), function(g) {
      context$trait <- M[g, ]
      call_fun(context, g)[1]
    }, numeric(1))
  }
  n_failed <- sum(is.na(values))
  if (n_failed) message("PPI ", entry$name, ": no value for ", n_failed, " genotype(s)")
  values
}

# PPI value of every genotype from its averaged reaction norm: the replicates are
# averaged per genotype and environment first (named numeric vector)
ppi_per_avg_genotype <- function(entry, phenotypes, environments) {
  compute_ppi(entry, replicate_means(phenotypes, environments), environments)
}

# PPI values of every genotype (rows) and replicate (columns) of one shape: the PPI is
# computed on each replicate's own measured reaction norm, nothing is averaged
ppi_per_replicate <- function(entry, phenotypes, environments) {
  ids <- unique(as.character(phenotypes$Genotype))
  reps <- sort(unique(phenotypes$Replicate))
  out <- sapply(reps, function(r) {
    d <- phenotypes[phenotypes$Replicate == r, ]
    if (anyDuplicated(paste(d$Genotype, d$Environment))) {
      stop("replicate ", r, ": more than one value per genotype and environment")
    }
    # the measured reaction norms of this replicate: one row per genotype, one column per environment
    M <- matrix(NA_real_, length(ids), length(environments), dimnames = list(ids, NULL))
    M[cbind(match(as.character(d$Genotype), ids), match(d$Environment, environments))] <- d$Trait
    compute_ppi(entry, M, environments)
  })
  matrix(out, length(ids), length(reps), dimnames = list(ids, reps))
}

# phenotypes: named list (one table per shape); environments: one vector for all
# shapes or a named list per shape (trait_environments()). Returns a named list (one table per
# PPI) with columns Score, Type (shape), Genotype_ID (line ID) and Replicate.
# replicate_handling decides on which level a PPI value is computed:
#   "mean_before"  replicates are averaged per genotype and environment first,
#                  then one PPI value per genotype (Replicate = NA); the standard
#   "replicates"   one PPI value per replicate's reaction norm, kept as separate rows
# A population-level PPI (FW) needs the reaction norms of the whole population and
# gives one value per genotype; it is always computed once on the replicate means,
# also with "replicates" (copies of the same value per replicate would only inflate n).
compute_ppis <- function(phenotypes, environments, registry = default_ppi_registry(), exclude = NULL,
                         replicate_handling = c("mean_before", "replicates")) {
  replicate_handling <- match.arg(replicate_handling)
  unknown <- setdiff(exclude, names(registry))
  if (length(unknown)) stop("cannot exclude unknown PPI(s): ", paste(unknown, collapse = ", "))
  registry <- registry[setdiff(names(registry), exclude)]
  lapply(registry, function(entry) {
    do.call(rbind, lapply(names(phenotypes), function(shape) {
      ph <- phenotypes[[shape]]
      env <- trait_environments(environments, shape)
      started <- Sys.time()
      on.exit(message(sprintf("PPI %s, %s: %d genotypes, %.1f s", entry$name, shape, length(unique(ph$Genotype)),
                              as.numeric(Sys.time() - started, units = "secs"))))
      if (replicate_handling == "mean_before" || entry$level == "population") {
        v <- ppi_per_avg_genotype(entry, ph, env)
        return(data.frame(Score = unname(v), Type = shape, Genotype_ID = names(v), Replicate = NA_integer_,
                          stringsAsFactors = FALSE))
      }
      P <- ppi_per_replicate(entry, ph, env)
      data.frame(Score = as.vector(t(P)), Type = shape, Genotype_ID = rep(rownames(P), each = ncol(P)),
                 Replicate = rep(as.integer(colnames(P)), times = nrow(P)), stringsAsFactors = FALSE)
    }))
  })
}

# PPI results (compute_ppis(), replicates averaged) as one table: one row per genotype
# (Shape, Genotype), one column per PPI; rows ordered by `shapes`, then as computed
ppi_wide_table <- function(traits, shapes = unique(traits[[1]]$Type)) {
  first <- traits[[1]]
  if (any(!is.na(first$Replicate))) stop("ppi_wide_table needs one value per genotype (replicate_handling = \"mean_before\")")
  key <- unique(first[, c("Type", "Genotype_ID")])
  key <- key[order(match(key$Type, shapes)), ]
  out <- data.frame(Shape = key$Type, Genotype = key$Genotype_ID, stringsAsFactors = FALSE)
  for (p in names(traits)) {
    d <- traits[[p]]
    out[[p]] <- d$Score[match(paste(out$Shape, out$Genotype), paste(d$Type, d$Genotype_ID))]
  }
  rownames(out) <- NULL
  out
}
