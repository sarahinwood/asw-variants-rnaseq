library(data.table)

snps_only <- fread("output/02_merged/missing_filtered_snps.out", header=F)
setnames(snps_only, old="V2", new="snpsonly_missing")
snps_only_no_variants <- 6771
snps_only$`%_snps_missing` <- (snps_only$snpsonly_missing/snps_only_no_variants)*100