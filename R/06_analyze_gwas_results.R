# =============================================================================
# 06_analyze_gwas_results.R — summarise simulated-GWAS results
# =============================================================================
# WHAT IT DOES: Loads the causal-SNP ground truth and the simulated GWAS result tables
#   and summarises, per plasticity index, how the recovered signal compares to truth,
#   rendering summary heatmaps/plots.
# REQUIRES:     GWAS result CSVs produced by 04/05 (run a scenario first); the causal
#               truth text files under the scenario's OUTPUT_BASE.
# PRODUCES:     Summary heatmaps/plots under OUTPUT_BASE.
# HOW TO RUN:   setwd("~/PlastQuant"); source(here::here("R", "06_analyze_gwas_results.R"))
# -----------------------------------------------------------------------------
# PARAMETERS (edit at the noted line)
#   PPIs   character vector   which plasticity indices to summarise; edit `PPIs <-`     [COMMON]
#   GWAS   logical            whether GWAS results are present
# =============================================================================
options(warn = -1)  # silence warnings even if the project .Rprofile was not loaded
## Import causal snp truth data from text files
library(pheatmap)

causalsnpstruth <- read.csv(
  "../causal_snps_truth.csv",
  header = FALSE,
  colClasses = c("character", "factor", "numeric"),
  col.names = c("numSNPS_SNP", "Parameter", "numeric"),
  stringsAsFactors = TRUE
)

causalsnpstruth$run = sub("SNP.*","",causalsnpstruth$numSNPS_SNP)
unique_run_params = unique(data.frame(causalsnpstruth$run, causalsnpstruth$Parameter))
for(i in seq_along(unique_run_params[,1])) {
  current_param_snps = (causalsnpstruth$run == unique_run_params$causalsnpstruth.run[i] & 
                          causalsnpstruth$Parameter == unique_run_params$causalsnpstruth.Parameter[i])
  causalsnpstruth$scaled[current_param_snps] = scale(causalsnpstruth$numeric[current_param_snps])
  causalsnpstruth$total[current_param_snps] = sum(abs(causalsnpstruth$scaled[current_param_snps]))
}

## ---- Import GWAS result data ---
gwasresults3plus <- read.csv(
  "../gwas_results2_significant.csv",
  header = FALSE,
  colClasses = c("numeric", "character", "numeric", "factor", "factor", "numeric", "numeric"),
  col.names = c("Heritability_pct", "numSNPS_SNP", "pVal", "PPI", "RN_type", "log10_pVar", "VarName7"),
  stringsAsFactors = TRUE
)

## ---- Define constants ----
RN_param_pairs <- matrix(c(
  "linear","Slope",
  "linear","BaseShift",
  "gaussian","Slope",
  "gaussian","BaseShift",
  "gaussian","Amplitude",
  "gaussian","Width",
  "gaussian","Center",
  "sinusoidal","Slope",
  "sinusoidal","BaseShift",
  "sinusoidal","Amplitude",
  "sinusoidal","Frequency",
  "sinusoidal","Phase",
  "wave","Slope",
  "wave","BaseShift",
  "wave","Amplitude",
  "wave","Frequency"
), ncol = 2, byrow = TRUE)

RN_type_indices = list(linear=c(1,2), gaussian=c(3:7), sinusoidal=c(8:12), wave=c(13:16))

PPIs <- c(
  "APC","CEV","CV_t","D_slope","ESP","ESPI","ESPIID",
  "EVS","gPi","NRW","PFI","PILSM","PImd","PIR",
  "PPF","PPi","PQ","PR","PSI","RC","RDPI",
  "RN","RNN","RPI","RSI","RTR","SI","FW"
)

numSNPs <- c("5","10","20")
H_vals <- c(20,40,60)

## ---- Initialize arrays ----
dims <- c(
  nrow(RN_param_pairs),
  length(PPIs),
  length(numSNPs),
  length(H_vals)
)

foundSNPsNum <- array(0, dim = dims)
foundCausalSNPsNum <- array(0, dim = dims)

foundSNPsLists <- vector("list", prod(dims))
dim(foundSNPsLists) <- dims

foundCausalSNPsLists <- vector("list", prod(dims))
dim(foundCausalSNPsLists) <- dims

## ---- Main analysis loop ----
for (PPI_i in seq_along(PPIs)) {
  for (nSNP_i in seq_along(numSNPs)) {
    for (H_i in seq_along(H_vals)) {
      for (rnp_i in seq_len(nrow(RN_param_pairs))) {
          ## Filter GWAS results
        subset_rows <- 
          gwasresults3plus[
            gwasresults3plus$PPI == PPIs[PPI_i] &
              gwasresults3plus$RN_type ==
              RN_param_pairs[rnp_i,1] &
              gwasresults3plus$Heritability_pct ==
              H_vals[H_i] &
              startsWith(
                gwasresults3plus$numSNPS_SNP,
                numSNPs[nSNP_i]
              ),
          ]
        
        foundSNPsLists[[rnp_i,PPI_i,nSNP_i,H_i]] <- subset_rows
        foundSNPsNum[rnp_i,PPI_i,nSNP_i,H_i] <- nrow(subset_rows)
        
        ## Identify causal SNPs
        current_causalSNPs <- causalsnpstruth$numSNPS_SNP[causalsnpstruth$Parameter == RN_param_pairs[rnp_i,2]]
        
        current_list <- subset_rows
        causal_rows <- current_list[current_list$numSNPS_SNP %in% current_causalSNPs, ]
        if (nrow(causal_rows) > 0) {
          foundCausalSNPsLists[[rnp_i,PPI_i,nSNP_i,H_i]] <- causal_rows
          foundCausalSNPsNum[rnp_i,PPI_i,nSNP_i,H_i] <- nrow(causal_rows)
        }
        

      }
      ## TODO store false positives already here!, but they need to be per Reaction_norm, not per reaction-norm-parameter!
      
    }
  }
}

foundCausalSNPs <- data.frame(matrix(NA, ncol=6, nrow = sum(foundCausalSNPsNum)))
colnames(foundCausalSNPs) = c("rnp_i", "PPI_i", "nSNP_i", "H_i", "numSNPS_SNP", "effect");
foundSNPs <- data.frame(matrix(NA, ncol=6, nrow = sum(foundSNPsNum)))
colnames(foundSNPs) = c("rnp_i", "PPI_i", "nSNP_i", "H_i", "numSNPS_SNP", "effect");
foundRN_SNPs <- data.frame(matrix(NA, ncol=6, nrow = sum(foundSNPsNum)))
colnames(foundRN_SNPs) = c("rnp_i", "PPI_i", "nSNP_i", "H_i", "numSNPS_SNP", "effect");


current_causal_entry = 0;
current_entry = 0;
current_RN_entry = 0;
for (PPI_i in seq_along(PPIs)) {
  for (nSNP_i in seq_along(numSNPs)) {
    for (H_i in seq_along(H_vals)) {
      for (rnp_i in seq_len(nrow(RN_param_pairs))) {
        if (foundCausalSNPsNum[rnp_i,PPI_i,nSNP_i,H_i] > 0) {
          for (row_i in seq (foundCausalSNPsNum[rnp_i,PPI_i,nSNP_i,H_i])) {
            current_causal_entry = current_causal_entry + 1;
            foundCausalSNPs$rnp_i[current_causal_entry] = rnp_i;
            foundCausalSNPs$PPI_i[current_causal_entry] = PPI_i;
            foundCausalSNPs$nSNP_i[current_causal_entry] = nSNP_i;
            foundCausalSNPs$H_i[current_causal_entry] = H_i;
            foundCausalSNPs$numSNPS_SNP[current_causal_entry] = foundCausalSNPsLists[[rnp_i, PPI_i, nSNP_i, H_i]]$numSNPS_SNP[row_i];
            foundCausalSNPs$effect[current_causal_entry] = 
              causalsnpstruth$scaled[causalsnpstruth$numSNPS_SNP == foundCausalSNPs$numSNPS_SNP[current_causal_entry]];
          }
        }
        if (foundSNPsNum[rnp_i,PPI_i,nSNP_i,H_i] > 0) {
          for (row_i in seq (foundSNPsNum[rnp_i,PPI_i,nSNP_i,H_i])) {
            current_entry = current_entry + 1;
            foundSNPs$rnp_i[current_entry] = rnp_i;
            foundSNPs$PPI_i[current_entry] = PPI_i;
            foundSNPs$nSNP_i[current_entry] = nSNP_i;
            foundSNPs$H_i[current_entry] = H_i;
            foundSNPs$numSNPS_SNP[current_entry] = foundSNPsLists[[rnp_i, PPI_i, nSNP_i, H_i]]$numSNPS_SNP[row_i];
            foundSNPs$effect[current_entry] = 0;
            if (sum(causalsnpstruth$numSNPS_SNP == foundSNPs$numSNPS_SNP[current_entry]) > 0) {
              # TODO correct formula (with scale (standardized per scenario))
              foundSNPs$effect[current_entry] = causalsnpstruth$scaled[causalsnpstruth$numSNPS_SNP == foundSNPs$numSNPS_SNP[current_entry]];
            }
            if (rnp_i %in% c(1,3,8,13)) {
              #foundSNPS are identical for all parameters of one simulated RN type
              #causal SNPs are only found for the exact RN type/parameter combination
              current_RN_entry = current_RN_entry + 1;
              foundRN_SNPs$rnp_i[current_RN_entry] = rnp_i;
              foundRN_SNPs$PPI_i[current_RN_entry] = PPI_i;
              foundRN_SNPs$nSNP_i[current_RN_entry] = nSNP_i;
              foundRN_SNPs$H_i[current_RN_entry] = H_i;
              foundRN_SNPs$numSNPS_SNP[current_RN_entry] = foundSNPsLists[[rnp_i, PPI_i, nSNP_i, H_i]]$numSNPS_SNP[row_i];
              foundRN_SNPs$effect[current_RN_entry] = 0;
              if (sum(causalsnpstruth$numSNPS_SNP == foundRN_SNPs$numSNPS_SNP[current_RN_entry]) > 0) {
                foundRN_SNPs$effect[current_RN_entry] = causalsnpstruth$scaled[causalsnpstruth$numSNPS_SNP == foundRN_SNPs$numSNPS_SNP[current_RN_entry]];
              }

            }
          }
        }
      }
    }
  }
}
foundRN_SNPs = foundRN_SNPs[1:current_RN_entry,]

foundCausalSNPsNoPPI = data.frame(foundCausalSNPs$numSNPS_SNP, foundCausalSNPs$H_i, foundCausalSNPs$rnp_i)
uniqueFoundCausalSNPs = unique(foundCausalSNPsNoPPI)

foundFalseSNPs = foundRN_SNPs[foundRN_SNPs$effect==0,]
foundFalseSNPsNoPPI = data.frame(foundFalseSNPs$numSNPS_SNP, foundFalseSNPs$H_i, foundFalseSNPs$rnp_i)
uniqueFoundFalseSNPs = unique(foundFalseSNPsNoPPI)

ppiFoundCausalSNPsMatrix = matrix(FALSE, ncol=length(PPIs), nrow = nrow(uniqueFoundCausalSNPs))
ppiFoundFalseSNPsMatrix = matrix(FALSE, ncol=length(PPIs), nrow = nrow(uniqueFoundFalseSNPs))

for (iCausal in 1:nrow(uniqueFoundCausalSNPs)) {
  ppiFoundCausalSNPsMatrix[iCausal, foundCausalSNPs$PPI_i[foundCausalSNPs$rnp_i == uniqueFoundCausalSNPs$foundCausalSNPs.rnp_i[iCausal] &  
                             foundCausalSNPs$H_i == uniqueFoundCausalSNPs$foundCausalSNPs.H_i[iCausal] & 
                             foundCausalSNPs$numSNPS_SNP == uniqueFoundCausalSNPs$foundCausalSNPs.numSNPS_SNP[iCausal]]] = TRUE
}

for (iFalse in 1:nrow(uniqueFoundFalseSNPs)) {
  ppiFoundFalseSNPsMatrix[iFalse, foundFalseSNPs$PPI_i[foundFalseSNPs$rnp_i == uniqueFoundFalseSNPs$foundFalseSNPs.rnp_i[iFalse] &  
                            foundFalseSNPs$H_i == uniqueFoundFalseSNPs$foundFalseSNPs.H_i[iFalse] & 
                            foundFalseSNPs$numSNPS_SNP == uniqueFoundFalseSNPs$foundFalseSNPs.numSNPS_SNP[iFalse]]] = TRUE
}

sum(rowSums(ppiFoundCausalSNPsMatrix[,c(22,26,28)])==0)
sum(rowSums(ppiFoundFalseSNPsMatrix[,c(22,26,28)])==0)

colSums(ppiFoundCausalSNPsMatrix[rowSums(ppiFoundCausalSNPsMatrix[,c(22,26,28)])==0,])
colSums(ppiFoundCausalSNPsMatrix[rowSums(ppiFoundCausalSNPsMatrix[,c(19,22,26,28)])==0,])

colSums(ppiFoundCausalSNPsMatrix[rowSums(ppiFoundCausalSNPsMatrix[,c(2,4,19,22,23,26,28)])==0,])



## ---- Aggregation ----

causalSNPSum_per_PPI <-
  apply(
    foundCausalSNPsNum,
    2,
    sum
  )

PPIs_sorted_i <-
  order(
    causalSNPSum_per_PPI,
    decreasing = TRUE
  )

causalSNPSum_per_RN_par <-
  apply(
    foundCausalSNPsNum,
    1,
    sum
  )

RN_par_sorted_i <-
  order(
    causalSNPSum_per_RN_par,
    decreasing = TRUE
  )


## ---- Create sorted 2D matrix ----

rows_out <-
  dim(foundCausalSNPsNum)[1] *
  dim(foundCausalSNPsNum)[3]

cols_out <-
  dim(foundCausalSNPsNum)[2] *
  dim(foundCausalSNPsNum)[4]

results_sorted_2d <-
  matrix(
    0,
    nrow = rows_out,
    ncol = cols_out
  )
results_unsorted_2d <-
  matrix(
    0,
    nrow = rows_out,
    ncol = cols_out
  )

results_causal_sorted_2d <-
  matrix(
    0,
    nrow = rows_out,
    ncol = cols_out
  )

results_causal_fractions_sorted_2d <-
  matrix(
    0,
    nrow = rows_out,
    ncol = cols_out
  )

results_causal_unsorted_2d <-
  matrix(
    0,
    nrow = rows_out,
    ncol = cols_out
  )

## ---- Labels ----

labels_Y <- character(rows_out)
labels_unsorted_Y <- character(rows_out)

for (rnp_i in seq_len(nrow(RN_param_pairs))) {
  for (nSNP_i in seq_along(numSNPs)) {
    idx <- 3*(rnp_i-1) + nSNP_i
    labels_Y[idx] <-
      paste(
        RN_param_pairs[
          RN_par_sorted_i[rnp_i],
          1
        ],
        RN_param_pairs[
          RN_par_sorted_i[rnp_i],
          2
        ],
        numSNPs[nSNP_i]
      )
    labels_unsorted_Y<-
      paste(
        RN_param_pairs[rnp_i,1],
        RN_param_pairs[rnp_i,2],
        numSNPs[nSNP_i]
      )
  }
}


labels_X <- character(cols_out)
labels_unsorted_X <- character(cols_out)

for (PPI_i in seq_along(PPIs)) {
  for (H_i in seq_along(H_vals)) {
    idx <- 3*(PPI_i-1) + H_i
    labels_X[idx] <- paste(PPIs[PPIs_sorted_i[PPI_i]], H_vals[H_i], "% H")
    labels_unsorted_X[idx] <- paste(PPIs[PPI_i], H_vals[H_i], "% H")
  }
}

## ---- Fill 2D matrix ----
for (rnp_i in seq_len(nrow(RN_param_pairs))) {
  for (PPI_i in seq_along(PPIs)) {
    for (nSNP_i in seq_along(numSNPs)) {
      for (H_i in seq_along(H_vals)) {
        row_idx <- 3*(rnp_i-1) + nSNP_i
        col_idx <- 3*(PPI_i-1) + H_i
        results_causal_sorted_2d[row_idx, col_idx] <-
          foundCausalSNPsNum[ RN_par_sorted_i[rnp_i], PPIs_sorted_i[PPI_i], nSNP_i, H_i]
        results_causal_fractions_sorted_2d[row_idx, col_idx] <-
          foundCausalSNPsNum[ RN_par_sorted_i[rnp_i], PPIs_sorted_i[PPI_i], nSNP_i, H_i]/(as.numeric(numSNPs[nSNP_i])*10)
        results_causal_unsorted_2d[row_idx, col_idx] <-
          foundCausalSNPsNum[ rnp_i, PPI_i, nSNP_i, H_i]
        results_sorted_2d[row_idx, col_idx] <-
          foundSNPsNum[ RN_par_sorted_i[rnp_i], PPIs_sorted_i[PPI_i], nSNP_i, H_i]
        results_unsorted_2d[row_idx, col_idx] <-
          foundSNPsNum[ rnp_i, PPI_i, nSNP_i, H_i]
      }
    }
  }
}


## ---- Plot heatmap ----

pheatmap(
  results_sorted_2d,
  labels_row = labels_Y,
  labels_col = labels_X,
  cluster_rows = FALSE,
  cluster_cols = FALSE,
)

pheatmap(
  results_causal_sorted_2d,
  labels_row = labels_Y,
  labels_col = labels_X,
  cluster_rows = FALSE,
  cluster_cols = FALSE,
)

pheatmap(
  results_causal_fractions_sorted_2d,
  labels_row = labels_Y,
  labels_col = labels_X,
  cluster_rows = FALSE,
  cluster_cols = FALSE,
)

pheatmap(
  results_unsorted_2d,
  labels_row = labels_Y,
  labels_col = labels_X,
  cluster_rows = FALSE,
  cluster_cols = FALSE,
)

pheatmap(
  results_causal_unsorted_2d,
  labels_row = labels_Y,
  labels_col = labels_X,
  cluster_rows = FALSE,
  cluster_cols = FALSE,
)
 
# title(
#   "causal SNPs with p<1E-3 per scenario, reaction norm parameter, and plasticity index"
# )


## ---- Histogram ----

hist(
  foundSNPsNum,
  main =
    "histogram of all SNPs with p<1E-3 per scenario, reaction norm parameter, and plasticity index"
)


hist(
  foundCausalSNPsNum,
  main =
    "histogram of all SNPs with p<1E-3 per scenario, reaction norm parameter, and plasticity index"
)

#hist(
#  foundCausalSNPsNum[foundCausalSNPsNum > 0],0.5:(max(foundCausalSNPsNum)+0.5),
#  main = "histogram of all SNPs with p<1E-3 per scenario, reaction norm parameter, and plasticity index"
#)

library(Matrix)

# suppose your sparse matrix is called m

# extract non-zero entries
sparseCausalResults <- Matrix(results_causal_sorted_2d, sparse = TRUE)
s = summary(sparseCausalResults)
sparseCausalData = data.frame(X <- s$j, Y <- s$i, Value <- s$x)   # values;

library(ggplot2)

ggplot(sparseCausalData, aes(x = X, y = Y, size = Value)) +
  geom_point(shape = 21, fill = "steelblue", alpha = 0.6) +
  scale_size_area(max_size = 10) +
  coord_fixed() +
  theme_minimal()

summedCausalData = rowSums(foundCausalSNPsNum, dims = 2);

s = summary(Matrix(summedCausalData, sparse = TRUE));
sparseSumData = data.frame(X <- s$j, Y <- s$i, Value <- s$x)   # values;
ggplot(sparseSumData, aes(x = X, y = Y, size = Value)) +
  geom_point(shape = 21, fill = "steelblue", alpha = 0.6) +
  scale_size_area(max_size = 10) +
  coord_fixed() +
  theme_minimal() +
  scale_x_continuous(label = PPIs, breaks=1:length(PPIs)) + 
  scale_y_continuous(label = apply(RN_param_pairs, 1, paste, collapse = " "), breaks=1:dim(RN_param_pairs)[1])

allPositives = foundSNPsNum;
allPositivesPerPPIScen = rowSums(allPositives, dims = 2)
allPositivesPerPPI = colSums(allPositivesPerPPIScen)
allPositivesPerScen = rowSums(allPositivesPerPPIScen)

allPositivesPerPPISlope = colSums(allPositivesPerPPIScen[c(1,3,8,13),])

allCausals = foundCausalSNPsNum;
allCausals[,,1,]=50;
allCausals[,,2,]=100;
allCausals[,,3,]=200;
allCausalsPerPPIScen = rowSums(allCausals, dims = 2)
allCausalsPerPPI = colSums(allCausalsPerPPIScen)
allCausalsPerScen = rowSums(allCausalsPerPPIScen)

allCausalsPerPPISlope = colSums(allCausalsPerPPIScen[c(1,3,8,13),])
allCausalsPerPPINoSlope = colSums(allCausalsPerPPIScen[c(2,4:7,9:12,14:16),])
allCausalsPerPPINoSlopeOrShift = colSums(allCausalsPerPPIScen[c(5:7,10:12,15:16),])

truePositives = foundCausalSNPsNum
truePositivesPerPPIScen = rowSums(truePositives, dims = 2)
truePositivesPerPPI = colSums(truePositivesPerPPIScen)
truePositivesPerScen = rowSums(truePositivesPerPPIScen)

truePositivesPerPPISlope = colSums(truePositivesPerPPIScen[c(1,3,8,13),])
truePositivesPerPPINoSlope = colSums(truePositivesPerPPIScen[c(2,4:7,9:12,14:16),])


falseNegatives = allCausals - truePositives

falsePositivesPerPPI_RN = allPositivesPerPPIScen[c(1,3,8,13),];
falsePositivesPerPPI_RN[1,]=allPositivesPerPPIScen[1,]-truePositivesPerPPIScen[1,]-truePositivesPerPPIScen[2,];
falsePositivesPerPPI_RN[2,]=allPositivesPerPPIScen[3,]-truePositivesPerPPIScen[3,]-truePositivesPerPPIScen[4,]-truePositivesPerPPIScen[5,]-truePositivesPerPPIScen[6,]-truePositivesPerPPIScen[7,];
falsePositivesPerPPI_RN[3,]=allPositivesPerPPIScen[8,]-truePositivesPerPPIScen[8,]-truePositivesPerPPIScen[9,]-truePositivesPerPPIScen[10,]-truePositivesPerPPIScen[11,]-truePositivesPerPPIScen[12,];
falsePositivesPerPPI_RN[4,]=allPositivesPerPPIScen[13,]-truePositivesPerPPIScen[13,]-truePositivesPerPPIScen[14,]-truePositivesPerPPIScen[15,]-truePositivesPerPPIScen[16,];

falsePositivesPerPPIScen = falsePositivesPerPPI_RN[c(1,1,2,2,2,2,2,3,3,3,3,3,4,4,4,4),]
falsePositivesPerPPI = colSums(falsePositivesPerPPI_RN)
falsePositivesPerScen = rowSums(falsePositivesPerPPIScen)
falsePositivesPerRN = rowSums(falsePositivesPerPPI_RN)

falsePositivesPerPPISlope = colSums(falsePositivesPerPPIScen[c(1,3,8,13),])
falsePositivesPerPPINoSlope = colSums(falsePositivesPerPPIScen[c(2,4:7,9:12,14:16),])


precisionPerPPIScen = truePositivesPerPPIScen / (truePositivesPerPPIScen + falsePositivesPerPPIScen)
recallPerPPIScen = truePositivesPerPPIScen / allCausalsPerPPIScen
precisionPerPPI = truePositivesPerPPI / (truePositivesPerPPI + falsePositivesPerPPI)
recallPerPPI = truePositivesPerPPI / allCausalsPerPPI

precisionPerPPISlope = truePositivesPerPPISlope / (truePositivesPerPPISlope + falsePositivesPerPPISlope)
recallPerPPISlope = truePositivesPerPPISlope / allCausalsPerPPISlope
precisionPerPPINoSlope = truePositivesPerPPINoSlope / (truePositivesPerPPINoSlope + falsePositivesPerPPINoSlope)
recallPerPPINoSlope = truePositivesPerPPINoSlope / allCausalsPerPPINoSlope

# Source - https://stackoverflow.com/a/64298954
# Posted by Duck
# Retrieved 2026-05-18, License - CC BY-SA 4.0

#Scaling factor
sf <- max(precisionPerPPI)/max(recallPerPPI)
precisionAndRecallPerPPI = c(precisionPerPPI,recallPerPPI*sf)
dataPPI = c(PPIs,PPIs)
precisionOrRecall = dataPPI
precisionOrRecall[seq(28)] = "precision"
precisionOrRecall[seq(29,56)] = "recall"
recallPerPPIScaled = recallPerPPI*sf

val_df = data.frame(dataPPI, precisionAndRecallPerPPI, precisionOrRecall)

ggplot(val_df, aes(x = dataPPI, y = precisionAndRecallPerPPI, fill = precisionOrRecall)) +
  geom_col(position = "dodge") + 
  scale_y_continuous(sec.axis = ~ ./sf)

sfSlope <- max(precisionPerPPISlope)/max(recallPerPPISlope)
precisionAndRecallPerPPISlope = c(precisionPerPPISlope,recallPerPPISlope*sfSlope)
recallPerPPISlopeScaled = recallPerPPISlope*sfSlope

valSlope_df = data.frame(dataPPI, precisionAndRecallPerPPISlope, precisionOrRecall)

ggplot(valSlope_df, aes(x = dataPPI, y = precisionAndRecallPerPPISlope, fill = precisionOrRecall)) +
  geom_col(position = "dodge") + 
  scale_y_continuous(sec.axis = ~ ./sfSlope)


sfNoSlope <- max(precisionPerPPINoSlope)/max(recallPerPPINoSlope)
precisionAndRecallPerPPINoSlope = c(precisionPerPPINoSlope,recallPerPPINoSlope*sfNoSlope)
recallPerPPINoSlopeScaled = recallPerPPINoSlope*sfNoSlope

valNoSlope_df = data.frame(dataPPI, precisionAndRecallPerPPINoSlope, precisionOrRecall)

ggplot(valNoSlope_df, aes(x = dataPPI, y = precisionAndRecallPerPPINoSlope, fill = precisionOrRecall)) +
  geom_col(position = "dodge") + 
  scale_y_continuous(sec.axis = ~ ./sfNoSlope)





pheatmap(
  recallPerPPIScen,
  labels_row = apply(RN_param_pairs, 1, paste, collapse = " "),
  labels_col = PPIs,
  cluster_rows = FALSE,
  cluster_cols = FALSE,
)

pheatmap(
  precisionPerPPIScen,
  labels_row = apply(RN_param_pairs, 1, paste, collapse = " "),
  labels_col = PPIs,
  cluster_rows = FALSE,
  cluster_cols = FALSE,
)

pheatmap(
  truePositivesPerPPIScen,
  labels_row = apply(RN_param_pairs, 1, paste, collapse = " "),
  labels_col = PPIs,
  cluster_rows = FALSE,
  cluster_cols = FALSE,
)

pheatmap(
  falsePositivesPerPPIScen,
  labels_row = apply(RN_param_pairs, 1, paste, collapse = " "),
  labels_col = PPIs,
  cluster_rows = FALSE,
  cluster_cols = FALSE,
)

temp1 = rbind(truePositivesPerPPIScen, falsePositivesPerPPI_RN)
true_and_falsePositivesPerPPIScen = temp1[c(1,2,17,3:7,18,8:12,19,13:16,20),]

labels_row1 = apply(RN_param_pairs, 1, paste, collapse = " ")
labels_row2 = paste(RN_param_pairs[c(1,3,8,13),1], "false positives")
labels_row_temp = c(labels_row1, labels_row2)
labels_rows = labels_row_temp[c(1,2,17,3:7,18,8:12,19,13:16,20)]

# move FW to alphabetically correct place
resorted_PPI_indices = c(seq(8), 28, seq(9,27))
col1=hcl.colors(26, "RdYlBu")
col2=hcl.colors(400, "RdYlBu")
log_true_and_falsePositivesPerPPIScen = log(true_and_falsePositivesPerPPIScen[,resorted_PPI_indices])
log_true_and_falsePositivesPerPPIScen[log_true_and_falsePositivesPerPPIScen<0]=NA
pheatmap(
  log_true_and_falsePositivesPerPPIScen,
  labels_row = labels_rows,
  labels_col = PPIs[resorted_PPI_indices],
  cluster_rows = FALSE,
  cluster_cols = FALSE,
  fontsize = 20,
  gaps_row = c(3,9,15),
  gaps_col=seq(27),
  angle_col = 90,
  legend_breaks = log(c(1,4,10,25,50,100,200)),
  legend_labels = c(1,4,10,25,50,100,200),
#  color = c(col1[26:14],col2[200:1])
)


causalsnpstruthNumSNPs = as.numeric(sub("_","",substr(causalsnpstruth$numSNPS_SNP,1,2)))
# TODO: change to correct formula
causalSNPsAbsoluteEffectTmp = causalsnpstruth$scaled / causalsnpstruth$total
# TODO: change to correct formula
causalSNPsAbsoluteEffect = c(causalSNPsAbsoluteEffectTmp*20, causalSNPsAbsoluteEffectTmp*40, causalSNPsAbsoluteEffectTmp*60)

foundCausalSNPsAbsoluteEffect = foundCausalSNPs$effect * foundCausalSNPs$H_i *20 / c(5,10,20)[foundCausalSNPs$nSNP_i]

hist_all_rel_causals = hist(causalsnpstruth$scaled,seq(-3,3,0.2))
hist_found_rel_causals = hist(unique(foundCausalSNPs$effect),seq(-3,3,0.2))
hist_rel_found_ratio = hist_found_rel_causals$counts / hist_all_rel_causals$counts
barplot(hist_rel_found_ratio)

foundCausalSNPsNoSlopeEffect = foundCausalSNPs$effect[!foundCausalSNPs$rnp_i %in% c(1,3,8,12)]

hist_allCausals = hist(causalSNPsAbsoluteEffect, plot = TRUE, breaks = seq(-30,30,2))
hist_foundCausals = hist(unique(foundCausalSNPsAbsoluteEffect), plot = TRUE, breaks = seq(-30,30,2))
hist_foundRatio = hist_foundCausals$counts / hist_allCausals$counts
barplot(hist_foundRatio)

foundCausalSNPsNoSlopeEffect = foundCausalSNPs$effect[!foundCausalSNPs$rnp_i %in% c(1,3,8,13)]
hist_all_rel_causals_no_slope = hist(causalsnpstruth$scaled[causalsnpstruth$Parameter != "Slope"],seq(-3,3,0.2))

hist_foundCausalsNoSlope = hist(unique(foundCausalSNPsNoSlopeEffect), plot = TRUE, breaks = seq(-3,3,0.2))
hist_foundRatioNoSlope = hist_foundCausalsNoSlope$counts / hist_all_rel_causals_no_slope$counts
barplot(hist_foundRatioNoSlope)



# consensus: SNPs found by 2+ or|and connected PPIs: precision and recall



# precision/recall per RN


# PPI efficiency: minimal set of PPIs finding all found SNPs


# add high-effect recall to precisionAndRecallPerPPI, precisionOrRecall

total_found_causals = sum(rowSums(ppiFoundCausalSNPsMatrix)>0)
message(paste("Number of all found causal SNPs:",as.character(total_found_causals)))

total_found_causals_slope = sum(rowSums(ppiFoundCausalSNPsMatrix[uniqueFoundCausalSNPs$foundCausalSNPs.rnp_i %in% c(1,3,8,13),])>0)
message(paste("Number of all found causal SNPs for slope:",as.character(total_found_causals_slope),as.character(100*total_found_causals_slope/total_found_causals),"%"))

total_found_causals_slope_FW_RN_RTR = sum(rowSums(ppiFoundCausalSNPsMatrix[uniqueFoundCausalSNPs$foundCausalSNPs.rnp_i %in% c(1,3,8,13),c(22,26,28)])>0)
message(paste("Number of all found causal SNPs for slope with FW, RN or RTR:",as.character(total_found_causals_slope_FW_RN_RTR),
              as.character(100*total_found_causals_slope_FW_RN_RTR/total_found_causals_slope),"%"))

total_found_causals_none_slope = sum(rowSums(ppiFoundCausalSNPsMatrix[!uniqueFoundCausalSNPs$foundCausalSNPs.rnp_i %in% c(1,3,8,13),])>0)
message(paste("Number of all found causal SNPs for other parameters than slope:",as.character(total_found_causals_none_slope)))

#colSums(ppiFoundCausalSNPsMatrix[rowSums(ppiFoundCausalSNPsMatrix[,c(1,2,4,14,15,19,20,22,26,28)])==0,])
PPIs[c(1,4,19,23)]
total_other_found_causals_no_slope = sum(rowSums(ppiFoundCausalSNPsMatrix[!uniqueFoundCausalSNPs$foundCausalSNPs.rnp_i %in% c(1,3,8,13),c(1,4,19,23)])>0)

total_other_found_causals = sum(rowSums(ppiFoundCausalSNPsMatrix[,c(1,4,19,23)])>0)
total_other_found_false = sum(rowSums(ppiFoundFalseSNPsMatrix[,c(1,4,19,23)])>0)

sum(rowSums(ppiFoundCausalSNPsMatrix[,c(22,26,28)])>2)
sum(rowSums(ppiFoundFalseSNPsMatrix[,c(22,26,28)])>2)

