library(rehh)
library(data.table)

###have to run each chromosome separately - have over 200,000 so it's too janky unless using to test subset

#tom's function
RunIhs <- function(contig_vcf, min_pct){
  message(contig_vcf)
  my_hh <- data2haplohh(hap_file = contig_vcf,
                        polarize_vcf = FALSE,
                        min_perc_geno.mrk = min_pct,
                        vcf_reader = "vcfR",
                        chr.name = "ASW_TRINITY_DN5834_c0_g1")
  my_scan <- scan_hh(my_hh,
                     polarized = FALSE,
                     phased = FALSE,
                     discard_integration_at_border = TRUE)
  #return(ihh2ihs(my_scan, freqbin = 1))
  return(data.table(my_scan))
}


Ruakura_files <- list.files("output/03_filtered_separate/", pattern="R.*.vcf.gz", full.names = TRUE)
names(Ruakura_files) <- gsub(".*/.*/(.*).vcf.gz", "\\1", Ruakura_files)

Dunedin_files <- list.files("output/03_filtered_separate/", pattern="[DI].*.vcf.gz", full.names = TRUE)
names(Dunedin_files) <- gsub(".*/.*/(.*).vcf.gz", "\\1", Dunedin_files)

##running rehh
dunedin_scans <- rbindlist(lapply(Dunedin_files, RunIhs, min_pct=0.9))
ruakura_scans <- rbindlist(lapply(Ruakura_files, RunIhs, min_pct=0.9))

# compare
xpehh <- data.table(ies2xpehh(scan_pop1 = data.frame(dunedin_scans),
                              scan_pop2 = data.frame(ruakura_scans),
                              popname1 = "Dunedin",
                              popname2 = "Ruakura",
                              p.adjust.method = "BH",
                              include_freq = TRUE))
