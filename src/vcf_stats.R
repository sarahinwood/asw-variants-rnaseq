library(data.table)

##what do any of these stats mean/how to interpret?

pi.all <- fread('output/04_stats/nucl_diversity_1kb.windowed.pi', header=T)
fst.all <- fread('output/04_stats/neg_vs_pos_fst_1kb.windowed.weir.fst')
taj.all <- fread('output/04_stats/TajD_1kb.Tajima.D', header=T)

# pi - diversity, higher = more
hist(pi.all$PI,br=200)
boxplot(pi.all$PI,ylab="diversity")

plot(pi.all$BIN_START,pi.all$PI,xlab="position",ylab="diversity")

# tajD
  # positive tajD = balancing selection (multiple alleles maintained), sudden population contraction
  # negative = selective sweep (fixation of one mutation), population expansion after bottleneck
hist(taj.all$TajimaD,br=200)
plot(taj.all$BIN_START,taj.all$TajimaD,xlab="position",ylab="Tajima's D")

# FST - fixation index (0-1)
  # High FST = a considerable degree of differentiation among populations
      # variation explained by population structure, no diversity shared between
  # Low FST = no differentiation between populations
hist(fst.all$WEIGHTED_FST,br=200)
plot(fst.all$BIN_START,fst.all$WEIGHTED_FST,xlab="position",ylab="Tajima's D")


# genes with high FST = different between locations
# genes with high TajD = many variants - not selected for
  # probably interested in low TajD in Ruakura ASW?
# low pi = high selection = low TajD?

# SNPs & annots
trinotate <- fread("data/trinotate_annots.csv")
trinotate$edited_ID <- paste("ASW", trinotate$`#gene_id`, sep="_")

genes_snps <- unique(pi.all$CHROM)

genes_SNPs_annots <- subset(trinotate, edited_ID %in% genes_snps)
annotated_gene_SNPs_bx <- subset(genes_SNPs_annots, !(genes_SNPs_annots$sprot_Top_BLASTX_hit=="")) 
annotated_gene_SNPs_bp <- subset(genes_SNPs_annots, !(genes_SNPs_annots$sprot_Top_BLASTP_hit=="")) 

  