##
# prepare data, add GO and KEGG terms 
# johaGL 2021
##
library("dplyr")
library("tidyverse")

setwd("~/BulkAnalysis_plusNetwork/")

odir <- "data/"

genesinfo <- read.table("data/genesinfo.csv",sep="\t",header=T)

natmiVizOut <- paste0("natmiOut/Young_D7/",  
   "Network_exp_0_spe_0_det_0.2_top_0_signal_lrc2p_weight_mean/")

## as bulk rnaseq was used, 'cluster' and celltype are the same!
lr_df <- read.table(paste0(natmiVizOut,"Edges.csv"),sep=",",header=T)

dplyr::sample_n(lr_df,10)
colnames(lr_df)

# > summary(lr_df$Ligand.average.expression.value)
# Min.  1st Qu.   Median     Mean  3rd Qu.     Max. 
# 0.034    1.800    9.231   91.413   39.655 6033.696 
# > summary(lr_df$Receptor.average.expression.value)
# Min.   1st Qu.    Median      Mean   3rd Qu.      Max. 
# 0.0339    1.9007    8.8838   29.9446   37.1832 2062.6146 
# > summary(lr_df$Edge.average.expression.weight)
# Min.   1st Qu.    Median      Mean   3rd Qu.      Max. 
# 0.0       9.7      67.8    3498.8     503.7 1735602.3 
# > summary(lr_df$Ligand.derived.specificity.of.average.expression.value)
# Min.   1st Qu.    Median      Mean   3rd Qu.      Max. 
# 0.0001041 0.0923763 0.2651131 0.3623238 0.5803305 1.0000000 
# > summary(lr_df$Receptor.derived.specificity.of.average.expression.value)
# Min.   1st Qu.    Median      Mean   3rd Qu.      Max. 
# 0.0002535 0.1070957 0.2661628 0.3387904 0.4822131 1.0000000 
# > summary(lr_df$Edge.average.expression.derived.specificity)
# Min.   1st Qu.    Median      Mean   3rd Qu.      Max. 
# 0.0000002 0.0111042 0.0524927 0.1239200 0.1588310 1.0000000 
###

# filter by L and R specificity index (keep >= fixed minimal specificity)
minspec = 0.15
filtered_lr <- lr_df %>% 
  filter(Ligand.derived.specificity.of.average.expression.value >= minspec) %>%
  filter(Receptor.derived.specificity.of.average.expression.value >= minspec)
dim(filtered_lr)

library("clusterProfiler")
# establishing GO and KEGG on specific age_day TPM matrices?:
#  PROBLEM: this will re-run repeatedly the same genes, 
#       so solution: look for genes that change A LOT, so for
#       the others we conclude their expression is 'constant' across time points 
# custom strategy: 
#       build a 3D matrix for each age:
#             - x,y,z will be  genes * samples * time
#               we have 4 time, but not for all cells!
#               calculate delta between z_2 and z_1
#                 then between z_3  z_2 and 
#                 then between z_4  z_3
#               stock all these deltas in a temporary vector 
#               and pick absolute max delta (with its sign)

# get enrichment by age :

changes = list()
itpm <- read.table(paste0("data/meanTPMat",age,day,".txt"),sep='\t',header=T)

## repare separator 

# DO NOT set ensemble ids as rownames ! 
# print(all(length(unique(itpm$ensembl)) == length(itpm$ensembl)))
# tpmids <- itpm$ensembl
# names(tpmids) <- genesinfo[match(tpmids,genesinfo$Geneid),]$symbol



geneList <-  ##  numeric vector: tau? , log(TPM+1) 
names(geneList) <- as.character(d[,1])  ##  named vector


####### ===================================================================
##### appendix : 
# 1.
### check differences between tau index and NATMI specificity indexes:
# example by hand tau specificity from Yanai:
thing <- c(0, 8, 0, 0, 0, 2, 0, 2, 0, 0, 0, 0)
num <- sum(1-(thing/max(thing)))
den <- length(thing)-1
num/den
# tau on real data:

# I pick gene ENSMUSG00000023224, which in table lr_df is :
irow = itpm[itpm$ensembl=='ENSMUSG00000023224',c(-1)]
tau_pick = (sum(1-(irow/max(irow))))/(length(irow)-1)
spec_indexNATMI <- lr_df %>% filter(Ligand.symbol=='Serping1') %>%
  select(Sending.cluster,
         Ligand.derived.specificity.of.average.expression.value) %>% unique()
spec_indexNATMI
## consider completing this script by filtering meanTPMAgeDay matrices 
## based on Yanai's tau index.

# 2. see the very highly expressed genes (colagen...)
hist(lr_df$Ligand.average.expression.value[
lr_df$Ligand.average.expression.value > 1000])
VERYHIGHEXP <- lr_df[
  lr_df$Ligand.average.expression.value > 1000,]

dim(VERYHIGHEXP)
View(VERYHIGHEXP)


