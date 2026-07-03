grep ',[3-6]' scenario_*/synthetic_data/scores_output/run_2026-03-13_*/gwas_results_mlm_80.csv | grep -v 'e-' > gwas_results_3plus.csv
cat gwas_results_3plus.csv | sed 's/scenario_//; s/_/,"/; s/:"//' > gwas_results2_3plus.csv
for dn in $(echo scenario_*); do cp $dn/synthetic_data/genetics/causal_snps_truth.csv causal_snps_truth_${dn}.csv; done
for nsnp in 5 10 20; do cat causal_snps_truth_scenario_20_${nsnp}.csv | sed "s/SNP_/${nsnp}SNP_/;"; done | awk -F',' '{print $2","$1","$3}' | grep -v '"SNP"'> causal_snps_truth.csv
