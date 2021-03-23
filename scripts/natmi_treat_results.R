##
# prepare data to yield 
# adjacency list ???
##
library("dplyr")
library("tidyverse")

setwd("~/BulkAnalysis_plusNetwork/")

odir <- "graphmatrices/"

natmiVizOut <- "natmiD7/Network_exp_0_spe_0_det_0.2_top_0_signal_lrc2p_weight_mean/"

## as bulk rnaseq was used, 'cluster' and celltype are the same!
lrmat <- read.table(paste0(natmiVizOut,"Edges.csv"),sep=",",header=T)

dplyr::sample_n(lrmat,10)
colnames(lrmat)


# compare values with those sent in input (meanTPM)
age="Young"
day="D7"
youngD7 <- read.table(paste0("data/meanTPM",age,day,".txt"), sep='\t',
                         header=T)


# check if coverage is ok there:
metadata <- readRDS("data/metadata.rds")
coverage <- read.table("data/COVmatrix.csv", sep="\t", header=T)
rownames(coverage) <- coverage$Gene.ID
coverage$Gene.ID <- NULL
covCOLNAMES <- metadata[match(colnames(coverage), metadata$sample),]$newname
coverage <- coverage[,!is.na(covCOLNAMES)] # drop out bad samples

genesZeroCoverage <- coverage %>% filter(rowSums(coverage) == 0)

badgenes <- rownames(youngD7)[rownames(youngD7) %in% rownames(genesZeroCoverage)]
goodgenes <- rownames(youngD7)[!rownames(youngD7) %in% rownames(genesZeroCoverage)]

#> max(as.numeric(unlist(coverage)))
#[1] 270906.8

# > summary(lrmat$Ligand.average.expression.value)
# Min.  1st Qu.   Median     Mean  3rd Qu.     Max. 
# 0.034    1.800    9.231   91.413   39.655 6033.696 
# > summary(lrmat$Receptor.average.expression.value)
# Min.   1st Qu.    Median      Mean   3rd Qu.      Max. 
# 0.0339    1.9007    8.8838   29.9446   37.1832 2062.6146 
# > summary(lrmat$Ligand.derived.specificity.of.average.expression.value)
# Min.   1st Qu.    Median      Mean   3rd Qu.      Max. 
# 0.0001041 0.0923763 0.2651131 0.3623238 0.5803305 1.0000000 
# > summary(lrmat$Receptor.derived.specificity.of.average.expression.value)
# Min.   1st Qu.    Median      Mean   3rd Qu.      Max. 
# 0.0002535 0.1070957 0.2661628 0.3387904 0.4822131 1.0000000 
# > summary(lrmat$Edge.average.expression.weight)
# Min.   1st Qu.    Median      Mean   3rd Qu.      Max. 
# 0.0       9.7      67.8    3498.8     503.7 1735602.3 
# > summary(lrmat$Edge.average.expression.derived.specificity)
# Min.   1st Qu.    Median      Mean   3rd Qu.      Max. 
# 0.0000002 0.0111042 0.0524927 0.1239200 0.1588310 1.0000000 

hist(lrmat$Ligand.average.expression.value[
  lrmat$Ligand.average.expression.value > 3000])
VERYHIGHLIGEXP <- lrmat[
  lrmat$Ligand.average.expression.value > 3000,]

dim(VERYHIGHLIGEXP)
