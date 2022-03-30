library(tidyverse)
library(data.table)
library(viridis)

# read in metadata
sample_table <- fread("data/sample_table.csv")
sample_info <- sample_table[,c(1,2,6,9,10)]
##attack status
attack_status_info <- fread("../asw-evasion-rnaseq/data/sample_table.csv")
sample_attack <- attack_status_info[,c(1,14)]
attacked_samples <- subset(sample_attack, Attacked=="Y")
##add to sample info
sample_info$Attack <- ifelse(sample_info$sample_name %in% attacked_samples$sample_name, "Yes", "Not observed")

# read in data
pca <- fread("output/04_plink/filtered_snps_plink_pca.eigenvec", header = FALSE)
eigenval <- scan("output/04_plink/filtered_snps_plink_pca.eigenval")

# sort out the pca data
# remove extra sample name column
pca <- pca[,-1]
# set names
names(pca)[1] <- "sample_name"
names(pca)[2:ncol(pca)] <- paste0("PC", 1:(ncol(pca)-1))

# merge with sample info
pca_info <- merge(pca, sample_info)

# first convert to percentage variance explained
pve <- data.frame(PC = 1:20, pve = eigenval/sum(eigenval)*100)

# make plot
ggplot(pve, aes(PC, pve)) + geom_bar(stat = "identity")+
  ylab("Percentage variance explained") + theme_light()

# calculate the cumulative sum of the percentage variance explained
cumsum(pve$pve)

# plot pca - PC1 = LOCATION
ggplot(pca_info, aes(PC1, PC2, shape=Location, colour=Attack))+
  geom_point(size = 3)+
  scale_colour_viridis(discrete=T)+ 
  coord_equal()+
  theme_light()+
  xlab(paste0("PC1 (", signif(pve$pve[1], 3), "%)"))+
  ylab(paste0("PC2 (", signif(pve$pve[2], 3), "%)"))

# PC1 populations
pc1_pos <- subset(pca_info, PC1>0)
fwrite(list(pc1_pos$sample_name), "output/04_plink/pc1_positive_samples.txt")
pc1_neg <- subset(pca_info, PC1<0)
fwrite(list(pc1_neg$sample_name), "output/04_plink/pc1_negative_samples.txt")

# test other PC combos - only PC1 creates two subgroups of samples
ggplot(pca_info, aes(PC2, PC3, shape=Location, colour=Attack))+
  geom_point(size = 3)+
  scale_colour_viridis(discrete=T)+ 
  coord_equal()+
  theme_light()+
  xlab(paste0("PC2 (", signif(pve$pve[2], 2), "%)"))+
  ylab(paste0("PC3 (", signif(pve$pve[3], 3), "%)"))

long_pca_info <- melt(pca_info)

# plot all pcas
ggplot(long_pca_info, aes(sample_name, value, colour=Location))+
  geom_point()+
  scale_colour_viridis(discrete=T)+ 
  theme_light()+
  facet_wrap(~variable)

# plot pca - EXPERIMENT
ggplot(pca_info, aes(PC1, PC2, colour=Experiment))+
  geom_point(size = 3)+
  scale_colour_viridis(discrete=T)+ 
  coord_equal()+
  theme_light()+
  xlab(paste0("PC1 (", signif(pve$pve[1], 3), "%)"))+
  ylab(paste0("PC2 (", signif(pve$pve[2], 3), "%)"))

##what is DAPC - discriminant analysis of principle components?
