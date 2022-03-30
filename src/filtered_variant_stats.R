library(data.table)

variant_stats <- fread("output/merged/filtered_variant_stats.txt", header=F, na.strings=".")

setnames(variant_stats, old=c("V1", "V2", "V3", "V4", "V5", "V6", "V7"),
         new=c("FS", "SOR", "MQRankSum", "ReadPosRankSum", "QD", "MQ", "DP"))

##remove rows with NA - removes 3846 rows
  #no NAs in DP or SOR though as made them fail
variant_stats_nona <- na.omit(variant_stats)

########
## FS ## trinity < 30, all below
########

plot(density(variant_stats_nona$FS))
sum(variant_stats_nona$QD>30)

########
## QD ## trinity > 2, all above
########

plot(density(variant_stats_nona$QD))
sum(variant_stats_nona$QD<2)

########
## DP ## me > 10, still 1022309 below 10, not on GATK page though - better to use QD **
########

sum(variant_stats_nona$DP<10)

## Depth - remove lower than 10 and above 2*mean
plot(density(variant_stats_nona$DP))

#########
## SOR ## GATK recc. < 3 ***
#########

plot(density(variant_stats_nona$SOR))
sum(variant_stats_nona$SOR<3)

#########
## MQ ## GATK recc. > 40, all are
#########

plot(density(variant_stats_nona$MQ))
sum(variant_stats_nona$MQ>40)/length(variant_stats_nona$MQ)

###############
## MQRankSum ## GATK recc. between -&+ 12.5 - all pass
###############

plot(density(variant_stats_nona$MQRankSum))

####################
## ReadPosRankSum ## GATK recc. between -&+ 8.0 - all pass
####################

plot(density(variant_stats_nona$ReadPosRankSum))

#### could still re-filter for depth>10 and SOR, all other filters pass