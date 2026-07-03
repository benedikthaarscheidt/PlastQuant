## Import causal snp truth data from text files
library(pheatmap)

## ---- Import GWAS result data ---
gwasresultsSignificant <- read.csv(
  "../gwas_results_maize_significant.csv",
  header = FALSE,
  colClasses = c("character", "numeric", "factor", "factor", "numeric", "numeric"),
  col.names = c("SNP", "pVal", "PPI", "RN_type", "log10_pVar", "eff_pVal"),
  stringsAsFactors = TRUE
)

gwasresults3plus <- read.csv(
  "../gwas_results_maize_3plus.csv",
  header = FALSE,
  colClasses = c("character", "numeric", "factor", "factor", "numeric", "numeric"),
  col.names = c("SNP", "pVal", "PPI", "RN_type", "log10_pVar", "eff_pVal"),
  stringsAsFactors = TRUE
)

gwasresults5plus = gwasresults3plus[gwasresults3plus$log10_pVar>5,]

PPIs = c("APC", "CEV","CV_t", "D_slope","ESP","ESPI", "ESPIID", "EVS","FW","gPi",
         "NRW",  "PFI", "PILSM","PImd", "PIR","PPF","PPi","PQ", "PR", "PSI",
         "RC", "RDPI", "RN", "RNN",  "RPI", "RSI", "RTR", "SI")
SNPs = unique(gwasresultsSignificant$SNP)
SNPs5plus = unique(gwasresults5plus$SNP)
RNs = c("leafArea", "WUE", "biomass") #unique(gwasresultsSignificant$RN_type)

foundSNPs = gwasresultsSignificant;
for(i in 1:length(PPIs))
{
  foundSNPs$PPI_i[foundSNPs$PPI == PPIs[i]] = i
}
for(i in 1:length(SNPs))
{
  foundSNPs$SNP_i[foundSNPs$SNP == SNPs[i]] = i
}
for(i in 1:length(RNs))
{
  foundSNPs$RN_i[foundSNPs$RN_type == RNs[i]] = i
}

ppiFoundSNPsMatrix = matrix(FALSE, ncol=length(PPIs), nrow = length(SNPs))
ppiFoundSNPsMatrixBiomass = matrix(FALSE, ncol=length(PPIs), nrow = length(SNPs))
ppiFoundSNPsMatrixLeafArea = matrix(FALSE, ncol=length(PPIs), nrow = length(SNPs))
ppiFoundSNPsMatrixWUE = matrix(FALSE, ncol=length(PPIs), nrow = length(SNPs))

for (iSNP in 1:length(SNPs)) {
  ppiFoundSNPsMatrix[iSNP, foundSNPs$PPI_i[foundSNPs$SNP_i == iSNP]] = TRUE
  ppiFoundSNPsMatrixBiomass[iSNP, foundSNPs$PPI_i[foundSNPs$SNP_i == iSNP & foundSNPs$RN_i == 1]] = TRUE
  ppiFoundSNPsMatrixLeafArea[iSNP, foundSNPs$PPI_i[foundSNPs$SNP_i == iSNP & foundSNPs$RN_i == 2]] = TRUE
  ppiFoundSNPsMatrixWUE[iSNP, foundSNPs$PPI_i[foundSNPs$SNP_i == iSNP & foundSNPs$RN_i == 3]] = TRUE
}


git_path = "~/CRC_1644_Z2_GWAS_simple/A-Datasets/"
load(file.path(git_path, "Genotyping_matrix_Negro_et_al_2019.Rdata"))
mat2 <- mat2[biomass_data_raw$Variety, ]
mat2reduced = mat2[,colnames(mat2) %in% SNPs]
MAF = colSums(mat2reduced>0)/244



scores_list <- readRDS("~/CRC_1644_Z2_GWAS_simple/scenario_maize/synthetic_data/scores_output/run_2026-06-19_11-04-38/scores_list.rds")

scores_list_arr3 = as.numeric(scores_list_arr2)
dim(scores_list_arr3)=c(732,28)
cor_scores = cor(scores_list_arr3)
scores_list_arr2 = scores_list_arr[,1,]
scores_list_arr = array(unlist(scores_list),dim=c(732,3,28))

