# =============================================================================
# kendall_rank_corr.R — Kendall rank correlations between indices
# =============================================================================
# WHAT IT DOES: Ranks the genotypes under every plasticity index and computes the
#   Kendall rank correlation of every index pair, for the population at its own
#   sampling resolution.
# REQUIRES:     POPULATION_DIR with its PPIs (compute_population_ppis()).
# PRODUCES:     Rank-correlation tables under POPULATION_DIR/analyses/rank_flops/.
# HOW TO RUN:   setwd("~/PlastQuant")
#               POPULATION_DIR <- here::here("output", "populations", "ppi_comparison")
#               source(here::here("R","analysis","kendall_rank_corr.R"))
# -----------------------------------------------------------------------------
# PARAMETERS (edit the named assignment near the top of the script)
#   POPULATION_DIR        path, REQUIRED                  population folder               [COMMON]
#   score_names           character vector                which indices to correlate
# =============================================================================
options(warn = -1)  # silence warnings even if the project .Rprofile was not loaded

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(readr)
  library(purrr)
})

score_names <- c("CV_t","RN","RNN","D_slope","RC","gPi","PPF","PPi","PImd","PILSM","RTR","PIR","RDPI","ESPI","ESPIID","PSI","RPI","PQ","PR","NRW","ESP","CEV","PFI","APC","SI","RSI","EVS")

coerce_to_df <- function(x) {
  if (is.null(x)) return(tibble::tibble())
  if (inherits(x, "data.frame")) return(tibble::as_tibble(x))
  if (is.atomic(x) || is.list(x)) {
    v <- tryCatch(unlist(x, use.names = FALSE), error = function(e) NULL)
    if (is.null(v)) return(tibble::tibble())
    return(tibble::tibble(score = suppressWarnings(as.numeric(v))))
  }
  tibble::tibble()
}

pick_col <- function(df, prefer, fallback = character(), numeric_fallback = FALSE) {
  if (is.null(df) || ncol(as.data.frame(df)) == 0) return(NA_character_)
  nm <- names(df); if (is.null(nm)) return(NA_character_)
  low <- tolower(nm)
  w <- which(low %in% tolower(prefer))
  if (length(w) > 0) return(nm[w[1]])
  if (length(fallback) > 0) {
    w <- which(low %in% tolower(fallback))
    if (length(w) > 0) return(nm[w[1]])
  }
  if (numeric_fallback) {
    num_cols <- which(vapply(df, is.numeric, logical(1)))
    if (length(num_cols) > 0) return(nm[tail(num_cols, 1)])
  }
  NA_character_
}

std_index_long <- function(obj, index_name) {
  df <- coerce_to_df(obj)
  if (nrow(df) == 0) return(tibble::tibble(index = character(), genotype_id = character(), score = double()))
  type_col  <- pick_col(df, c("type","form","shape","genotype_form","gtype_form"))
  geno_col  <- pick_col(df, c("genotype","genotype_id","geno","id","name","sample","accession","line","individual"))
  score_col <- pick_col(df, c("score","value","index","metric"), numeric_fallback = TRUE)
  if (is.na(score_col)) return(tibble::tibble(index = character(), genotype_id = character(), score = double()))
  geno_vec <- if (!is.na(geno_col)) as.character(df[[geno_col]]) else as.character(seq_len(nrow(df)))
  gid <- if (!is.na(type_col)) paste0(tolower(as.character(df[[type_col]])),"::",geno_vec) else geno_vec
  sc <- suppressWarnings(as.numeric(df[[score_col]]))
  dplyr::tibble(index = index_name, genotype_id = gid, score = sc) |>
    dplyr::filter(!is.na(score), genotype_id != "")
}

pair_tau <- function(df_long, i, j) {
  a <- df_long |> dplyr::filter(index == i) |> dplyr::select(genotype_id, score) |> dplyr::distinct()
  b <- df_long |> dplyr::filter(index == j) |> dplyr::select(genotype_id, score) |> dplyr::distinct()
  ab <- dplyr::inner_join(a, b, by = "genotype_id", suffix = c("_i","_j")) |>
    tidyr::drop_na(score_i, score_j)
  est <- NA_real_; pval <- NA_real_
  if (nrow(ab) >= 3 && stats::sd(ab$score_i, na.rm = TRUE) > 0 && stats::sd(ab$score_j, na.rm = TRUE) > 0) {
    ct <- try(stats::cor.test(ab$score_i, ab$score_j, method = "kendall", exact = FALSE), silent = TRUE)
    if (!inherits(ct, "try-error")) { est <- unname(ct$estimate); pval <- unname(ct$p.value) }
  }
  dplyr::tibble(score_i = i, score_j = j, correlation = est, p_value = pval, n_common = nrow(ab))
}

if (!exists("POPULATION_DIR")) stop("set POPULATION_DIR (a population folder) first")
source(here::here("R", "pipeline", "load_pipeline.R"))
out_dir <- file.path(POPULATION_DIR, "analyses", "rank_flops")
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

# Scores of the population at its own resolution (one row per genotype)
score_tab <- population_score_table(POPULATION_DIR)
n_points <- nrow(read.csv(file.path(POPULATION_DIR, "data", "environments.csv")))
genotype_ids <- paste(score_tab$Shape, score_tab$Genotype, sep = "_")
scores_list <- list()
for (nm in intersect(score_names, names(score_tab))) {
  scores_list[[nm]] <- data.frame(genotype_id = genotype_ids, score = score_tab[[nm]])
}
if (length(scores_list) < 2) stop("Fewer than 2 scores in ", POPULATION_DIR)
x_long <- purrr::imap_dfr(scores_list, ~std_index_long(.x, .y))
if (nrow(x_long) == 0) stop("No valid score rows produced")
rank_mat <- x_long |>
  dplyr::group_by(index) |>
  dplyr::mutate(rank = base::rank(-score, ties.method = "average")) |>
  dplyr::ungroup() |>
  dplyr::select(genotype_id, index, rank) |>
  dplyr::distinct() |>
  tidyr::pivot_wider(names_from = index, values_from = rank) |>
  dplyr::arrange(genotype_id) |>
  dplyr::mutate(gid = dplyr::row_number()) |>
  dplyr::relocate(gid, .before = genotype_id)
idxs <- x_long |> dplyr::distinct(index) |> dplyr::arrange(index) |> dplyr::pull(index)
pairs <- if (length(idxs) >= 2) utils::combn(idxs, 2, simplify = FALSE) else list()
pairs_tbl <- purrr::map_dfr(pairs, ~pair_tau(x_long, .x[1], .x[2])) |>
  dplyr::mutate(resolution = n_points, .before = 1)
tag <- paste0("population_res", n_points)
readr::write_csv(pairs_tbl, file.path(out_dir, paste0("kendall_pairs_", tag, ".csv")))
readr::write_csv(rank_mat,  file.path(out_dir, paste0("ranks_matrix_combined_", tag, ".csv")))
readr::write_csv(rank_mat |> dplyr::select(gid, genotype_id), file.path(out_dir, paste0("genotype_map_", tag, ".csv")))
message("Saved: ", file.path(out_dir, paste0("kendall_pairs_", tag, ".csv")))
message("Saved: ", file.path(out_dir, paste0("ranks_matrix_combined_", tag, ".csv")))
message("Saved: ", file.path(out_dir, paste0("genotype_map_", tag, ".csv")))
