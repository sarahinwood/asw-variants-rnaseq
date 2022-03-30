library(data.table)
library(pcadapt)

##variant calling should be one sample at a time - not all

##skip lines with detail of headers and gene names
vcf <- fread("output/variants/IALM_NC_4/renamed_IALM_NC_4.vcf", skip = 239477)
trinotate <- fread("data/trinotate_annots.csv", na.strings = "")

##filter to keep only variants that passed filtering
filter_vcf <- subset(vcf, FILTER=="PASS")
##number genes with variants - 104,541
length(unique(filter_vcf$`#CHROM`))

##filter by read depth also
filter_vcf$DP <- tstrsplit(filter_vcf$INFO, "DP=", keep=c(2))
filter_vcf$DP <- tstrsplit(filter_vcf$DP, ";", keep=c(1))
filter_vcf$DP <- as.numeric(filter_vcf$DP)
##filter for read depth above 10 (Trinity/GATK suggestion)
depth_filter_vcf <- subset(filter_vcf, DP>10)
##number genes with variants - 85,020
length(unique(depth_filter_vcf$`#CHROM`))

##how many variants have trinity annotations
variant_annots <- subset(trinotate, `#gene_id` %in% depth_filter_vcf$`#CHROM`)
##13,140 have BlastX annotation
sum(!is.na(variant_annots$sprot_Top_BLASTX_hit))

##how to link back to samples?
