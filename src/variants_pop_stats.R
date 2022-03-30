library(data.table)
library(vcfR)
library(poppr)
library(hierfstat)
library(ggplot2)

metadata <- read.table("data/full_sample_table.csv", sep = ",", h=T)
metadata$Location <- as.factor(metadata$Location)

dataVCF <- read.vcfR(("output/02_merged/filtered_snps.vcf.gz"))
# convert to genind
geninddata<-vcfR2genind(dataVCF)
# assign populations
geninddata@pop<-(metadata[,2])

# basic stats
pop <- genind2hierfstat(geninddata, pop=metadata[,2])
basic.stats(pop)$overall
##can this be done on each population?

# pairwise fst - location only
pop_matrix_fst<-as.matrix(genet.dist(pop,method="WC84"))
colnames(pop_matrix_fst)<-levels(pop$pop)
rownames(pop_matrix_fst)<-levels(pop$pop)
pop_matrix_fst
pop_fst_long <- melt(pop_matrix_fst)
# plot pairwise fst location
ggplot(pop_fst_long, aes(x = Var1, y = Var2, fill = value)) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)) +
  xlab(NULL) + ylab(NULL) +
  scale_x_discrete(expand = c(0, 0)) +
  scale_y_discrete(expand = c(0, 0)) +
  scale_fill_viridis_c(option="A", direction=c(-1), limits=c(-0.0004,0.05), guide = guide_colorbar(title = expression(italic(F)["ST"])),
                       na.value = NA) +
  coord_fixed () +
  geom_raster()



# pairwise fst - location & attack
pop_attack <- genind2hierfstat(geninddata, pop=paste(metadata[,2],metadata[,6],sep="_"))
matrix_fst<-as.matrix(genet.dist(pop_attack,method="WC84"))
pop_attack$pop <- as.factor(pop_attack$pop)
colnames(matrix_fst)<-levels(pop_attack$pop)
rownames(matrix_fst)<-levels(pop_attack$pop)
matrix_fst

# make long & merge with groups
fst_long <- melt(matrix_fst)
sample_vs_pop <- metadata[,c(1,2,6)]
sample_vs_pop$group <- paste(sample_vs_pop$Location, sample_vs_pop$Attack_status, sep="_")
sample_group <- sample_vs_pop[,c(1,4)]
var1_long_fst <- merge(sample_group, fst_long, by.x=)
# plot fst loc_attack
ggplot(fst_long, aes(x = Var1, y = Var2, fill = value)) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)) +
  xlab(NULL) + ylab(NULL) +
  scale_x_discrete(expand = c(0, 0)) +
  scale_y_discrete(expand = c(0, 0)) +
  scale_fill_viridis_c(option="A", direction=c(-1), limits=c(-0.0004,0.05), guide = guide_colorbar(title = expression(italic(F)["ST"])),
    na.value = NA) +
  coord_fixed () +
  geom_raster()



# pairwise fst - all samples
pop_attack <- genind2hierfstat(geninddata, pop=paste(metadata[,1]))
matrix_fst<-as.matrix(genet.dist(pop_attack,method="WC84"))
pop_attack$pop <- as.factor(pop_attack$pop)
colnames(matrix_fst)<-levels(pop_attack$pop)
rownames(matrix_fst)<-levels(pop_attack$pop)
matrix_fst
