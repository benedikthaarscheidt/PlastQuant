# =============================================================================
# 01_simulate_genetics.R — simulate SNP genotypes and genetic parameters
# =============================================================================
# WHAT IT DOES: Generates synthetic SNP genotypes for NUM_GENOTYPES individuals with
#   optional nested population structure, assigns causal SNPs and effect sizes to the
#   reaction-norm parameters, and records the ground-truth causal map later used to
#   score GWAS accuracy.
# REQUIRES:     Nothing upstream (pure simulation; reads no external data).
# PRODUCES:     In memory: the genotype matrix and truth table consumed by 02/03.
#               If OUTPUT_BASE is set, writes them under OUTPUT_BASE.
# HOW TO RUN:   setwd("~/PlastQuant"); source(here::here("R", "01_simulate_genetics.R"))
# -----------------------------------------------------------------------------
# PARAMETERS (edit here)
#   NUM_GENOTYPES         integer, default 80      number of simulated genotypes           [COMMON]
#   GENETIC_VARIANCES     logical/vector, default FALSE  per-parameter genetic variances
#                                                  (FALSE = use built-in defaults)
#   POLY_MODE             "uniform"|"exponential", default "uniform"  causal-effect distribution
#   POLY_STRENGTH         integer 1-5, default 5   genetic-variance strength (5 = full target SD)
#   POLY_COUNTS           integer/vector, default 1  number of causal SNPs per parameter
#   STRUCTURED_POPULATION logical, default FALSE   add nested population structure
#   SEED                  integer, default 42      RNG seed for reproducibility            [COMMON]
#   OUTPUT_BASE           path, default here::here("output","default")  where outputs go
# -----------------------------------------------------------------------------
# These parameters are supplied by the sourcing script (03_plasticity_scores.R or a
# scenario_*.R driver) or by your R session before sourcing. They are intentionally
# NOT assigned here: the simulation branches on whether SEED / GENETIC_VARIANCES exist
# (see below), so defining them here would change the RNG path. GENETIC_VARIANCES falls
# back to FALSE further down; SEED, when defined, offsets the seeds (leave it unset to
# reproduce the default simulation).
# =============================================================================


pick_snps <- function(n, exclude, total_snps_num) {
  candidates <- setdiff(1:total_snps_num, exclude)
  if (length(candidates) < n) stop("Not enough SNPs left!")
  sample(candidates, n)
}

# Helper to calculate phenotype from SNPs
calc_pheno <- function(n_snps, mean_val, sd_val, total_snps_num, used_snps, snp_matrix, strength = 5, min_val = -Inf, max_val = Inf) {
  idx <- pick_snps(n_snps, used_snps, total_snps_num)

  # Generate Effects based on Mode
  if (exists("POLY_MODE") && POLY_MODE == "exponential") {
    # Exponential distribution (few strong, many weak)
    # rexp(n) gives positive values. Multiply by random sign.
    eff <- rexp(n_snps, rate = 1) * sample(c(-1, 1), n_snps, replace = TRUE)
  } else {
    # Uniform distribution (standard/box-like)
    eff <- runif(n_snps, -1, 1)
  }
  
  gv <- as.vector(snp_matrix[, idx, drop = FALSE] %*% eff)
  
  # Scale Genetic Variance based on Strength (1-5)
  # Strength 5 = 100% of target SD
  # Strength 1 = 20% of target SD
  effective_sd <- sd_val * (strength / 5)
  
  # Robust scaling
  if (sd(gv) == 0) {
    scaled <- rep(mean_val, length(gv))
  } else {
    scaled <- as.vector(scale(gv)) * effective_sd + mean_val
  }
  
  # Enforce bounds
  scaled[scaled < min_val] <- min_val
  scaled[scaled > max_val] <- max_val
  
  return(list(values = scaled, snps = idx, effects = eff, used_snps = used_snps))
}

create_truth_rows <- function(param_name, p_obj, snp_matrix) {
  # Check if the SNPs field is NULL or explicitly NA (i.e., no causal SNPs)
  if (is.null(p_obj$snps) || (length(p_obj$snps) == 1 && is.na(p_obj$snps[1]))) {
    # If no causal SNPs, return an empty data frame with the correct column structure.
    return(data.frame(Parameter = character(0), SNP = character(0), Effect = numeric(0)))
  }
  
  # If there are causal SNPs, return the data frame as usual
  data.frame(
    Parameter = param_name,
    SNP = colnames(snp_matrix)[p_obj$snps],
    Effect = p_obj$effects
  )
}

simulate_genetics <- function(num_forms = 4) {

if (exists("SEED")) {
  seed_add_on = SEED;
} else {
  seed_add_on = 0;
}

set.seed(123 + seed_add_on)

# -----------------------------------------------------------------------------
# 1. Configuration
# -----------------------------------------------------------------------------
if (exists("NUM_GENOTYPES")) {
  num_genotypes <- NUM_GENOTYPES
} else {
  num_genotypes <- 200
}

if (exists("NUM_SNPs")) {
  num_snps <- NUM_SNPs
} else {
  num_snps <- 10000
}

output_dir <- paste(OUTPUT_BASE, "/synthetic_data/genetics",sep="")

# Structure: 4 Forms. Dynamically calculate families and genotypes per family.

if (num_genotypes %% num_forms != 0) {
  stop("num_genotypes must be divisible by 4 (num_forms)")
}
genotypes_per_form <- num_genotypes / num_forms

# Try to keep genotypes_per_family around 10
if (genotypes_per_form %% 10 == 0) {
  genotypes_per_family <- 10
} else if (genotypes_per_form %% 5 == 0) {
  genotypes_per_family <- 5
} else {
  # Fallback: 1 family per form, all genotypes in it
  genotypes_per_family <- genotypes_per_form
}

families_per_form <- genotypes_per_form / genotypes_per_family
num_families <- num_forms * families_per_form
if (num_families * genotypes_per_family != num_genotypes) {
  stop("Configuration Error: num_families * genotypes_per_family must equal num_genotypes")
}
if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE)
}

# -----------------------------------------------------------------------------
# 2. Generate SNP Matrix with Population Structure
# -----------------------------------------------------------------------------
# We use a Balding-Nichols model:
# 1. Ancestral Allele Frequencies (P_anc) ~ Uniform(0.05, 0.5)
# 2. Family Allele Frequencies (P_fam) ~ Beta distribution parameterized by Fst
# 3. Genotypes ~ Binomial(2, P_fam)

# Set to FALSE for random mating (No Kinship needed)
Fst <- 0.2 # Only used if STRUCTURED_POPULATION is TRUE
# If Fst = 0: All families have the same allele frequencies (Random Mating).
# If Fst = 1: Families are completely fixed for different alleles (Different Species).
P_anc <- runif(num_snps, 0.2, 0.8)
# The baseline frequency of the minor allele for each SNP in the original "ancestral" population.
# Function to generate family frequencies
rbeta_fst <- function(n, p, fst) {
  alpha <- p * (1 - fst) / fst
  beta <- (1 - p) * (1 - fst) / fst
  return(rbeta(n, alpha, beta))
}

snp_matrix <- matrix(NA, nrow = num_genotypes, ncol = num_snps)
rownames(snp_matrix) <- paste0("Genotype_", 1:num_genotypes)
colnames(snp_matrix) <- paste0("SNP_", 1:num_snps)
#snp_matrix_per_RN = data.frame()
#for ( rn_form in c("linear", "gaussian", "sinusoidal", "wave")) {
#  tmp_matrix = matrix(NA, nrow = genotypes_per_form, ncol = num_snps)
#  rownames(tmp_matrix) <- paste0("Genotype_", 1:genotypes_per_form)
#  colnames(tmp_matrix) <- paste0("SNP_", 1:num_snps)
#  snp_matrix_per_RN[seq(genotypes_per_form), rn_form] = tmp_matrix
#}


if (STRUCTURED_POPULATION) {
  message("Generating Structured Population (Balding-Nichols Model)")
  # Generate data family by family
  for (f in 1:num_families) {
    # Generate allele frequencies for this family
    P_fam <- sapply(P_anc, function(p) rbeta_fst(1, p, Fst))

    # Indices for this family
    start_idx <- (f - 1) * genotypes_per_family + 1
    end_idx <- f * genotypes_per_family

    # Generate genotypes
    for (i in start_idx:end_idx) {
      snp_matrix[i, ] <- rbinom(num_snps, 2, P_fam)
    }
  }
} else {
  message("Generating Unstructured Population (Random Mating)")
  # Generate genotypes directly from ancestral frequencies
  for (i in 1:num_genotypes) {
    snp_matrix[i, ] <- rbinom(num_snps, 2, P_anc)
  }
}

# -----------------------------------------------------------------------------
# 3. Define Causal SNPs (Comprehensive Model)
# -----------------------------------------------------------------------------
used_snps <- c()


# --- Universal Parameters ---
# Slope: Mean 0, SD 3.0
p_slope <- calc_pheno(POLY_COUNTS$Slope, 0, 3.0, num_snps, used_snps, snp_matrix, strength = POLY_STRENGTH$Slope)
used_snps = c(used_snps, p_slope$idx)

# BaseShift: Mean 20, SD 5
p_shift <- calc_pheno(POLY_COUNTS$BaseShift, 20, 5, num_snps, used_snps, snp_matrix, strength = POLY_STRENGTH$BaseShift)
used_snps = c(used_snps, p_shift$idx)

# Check GENETIC_VARIANCES toggle (default to TRUE if not set, but user wants FALSE baseline)
if (!exists("GENETIC_VARIANCES")) GENETIC_VARIANCES <- FALSE

if (GENETIC_VARIANCES) {
  message("Modeling VE and VS as Genetic Traits (vQTLs)...")
  # VE (Noise): Mean 0.5, SD 0.1
  p_ve <- calc_pheno(POLY_COUNTS$VE, 0.5, 0.1, num_snps, used_snps, snp_matrix, min_val = 0.1, strength = POLY_STRENGTH$VE)
  used_snps = c(used_snps, p_ve$idx)
  
  # VS (Slope Var): Mean 0.2, SD 0.05
  p_vs <- calc_pheno(POLY_COUNTS$VS, 0.2, 0.05, num_snps, used_snps, snp_matrix, min_val = 0.01, strength = POLY_STRENGTH$VS)
  used_snps = c(used_snps, p_vs$idx)
} else {
  message("Modeling VE and VS as Environmental/Random (Non-Genetic)...")
  # Generate random values (not linked to SNPs)
  
  ve_vals <- rnorm(num_genotypes, 0.05, 0.01) # <<< MODIFIED FOR LOW NOISE
  ve_vals[ve_vals < 0.01] <- 0.01 
  p_ve <- list(values = ve_vals, snps = NA, effects = NA)
  
  # V_S (Slope Var): Target Mean 0.01, SD 0.005
  vs_vals <- rnorm(num_genotypes, 0.05, 0.01) # <<< MODIFIED FOR LOW NOISE
  vs_vals[vs_vals < 0.01] <- 0.01 # Ensure non-negative and minimum
  p_vs <- list(values = vs_vals, snps = NA, effects = NA)
}

# --- Shape-Specific Parameters ---
# Amplitude: Mean 5, SD 3, Min 1
p_amp <- calc_pheno(POLY_COUNTS$Amplitude, 5, 3, num_snps, used_snps, snp_matrix, min_val = 1, strength = POLY_STRENGTH$Amplitude)
used_snps = c(used_snps, p_amp$idx)

# Width: Mean 0.2, SD 0.05, Min 0.05
p_width <- calc_pheno(POLY_COUNTS$Width, 0.2, 0.05, num_snps, used_snps, snp_matrix, min_val = 0.05, strength = POLY_STRENGTH$Width)
used_snps = c(used_snps, p_width$idx)

# Center: Mean 5, SD 2
p_center <- calc_pheno(POLY_COUNTS$Center, 5, 2, num_snps, used_snps, snp_matrix, strength = POLY_STRENGTH$Center)
used_snps = c(used_snps, p_center$idx)

# Frequency: Mean 0.5, SD 0.1, Min 0.1
p_freq <- calc_pheno(POLY_COUNTS$Frequency, 0.5, 0.1, num_snps, used_snps, snp_matrix, min_val = 0.1, strength = POLY_STRENGTH$Frequency)
used_snps = c(used_snps, p_freq$idx)

# Phase: Mean 0, SD 1
p_phase <- calc_pheno(POLY_COUNTS$Phase, 0, 1, num_snps, used_snps, snp_matrix, strength = POLY_STRENGTH$Phase)
used_snps = c(used_snps, p_phase$idx)


# -----------------------------------------------------------------------------
# 4. Save Outputs
# -----------------------------------------------------------------------------

write.csv(snp_matrix, file.path(output_dir, "snp_matrix.csv"), row.names = TRUE)

params_df <- data.frame(
  Genotype_ID = 1:num_genotypes,
  Genotype_Name = rownames(snp_matrix),
  Form = rep(c("Linear", "Gaussian", "Sinusoidal", "Wave"), each = num_genotypes / 4),

  # Universal
  Slope = p_slope$values,
  BaseShift = p_shift$values,
  VE = p_ve$values,
  VS = p_vs$values,

  # Shape Specific
  Amplitude = p_amp$values,
  Width = p_width$values,
  Center = p_center$values,
  Frequency = p_freq$values,
  Phase = p_phase$values
)
write.csv(params_df, file.path(output_dir, "genotypic_parameters.csv"), row.names = FALSE)

# Truth Table

truth_df <- rbind(
  create_truth_rows("Slope", p_slope, snp_matrix),
  create_truth_rows("BaseShift", p_shift, snp_matrix),
  create_truth_rows("VE", p_ve, snp_matrix),
  create_truth_rows("VS", p_vs, snp_matrix),
  create_truth_rows("Amplitude", p_amp, snp_matrix),
  create_truth_rows("Width", p_width, snp_matrix),
  create_truth_rows("Center", p_center, snp_matrix),
  create_truth_rows("Frequency", p_freq, snp_matrix),
  create_truth_rows("Phase", p_phase, snp_matrix)
)
write.csv(truth_df, file.path(output_dir, "causal_snps_truth.csv"), row.names = FALSE)

cat("Simulation complete (80 Genotypes, Comprehensive Parameters).\n")
cat("Parameters saved to:", file.path(output_dir, "genotypic_parameters.csv"), "\n")

return(list(snp_matrix=snp_matrix, params_df=params_df, truth_df = truth_df))
}
