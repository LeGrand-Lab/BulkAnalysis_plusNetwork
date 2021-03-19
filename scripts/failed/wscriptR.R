library(dplyr)
library(tidyverse)
library(ggplot2)
library(cowplot)
library(forcats)
library(RColorBrewer)
library(MASS)
library(pheatmap)
library(DESeq2)
library(ggsci) # publishing palettes
library(cowplot)
library(gridExtra)
# devtools::install_github("vqv/ggbipplot")
# install.packages("glmpca")
# ** see
# latest : http://www.bioconductor.org/packages/release/bioc/vignettes/DESeq2/inst/doc/DESeq2.html#if-i-have-multiple-groups-should-i-run-all-together-or-split-into-pairs-of-groups
# **
### >>> there exist important technical variability between replicates !

setwd("~/bulk_analysis/")
fmat <- readRDS(prefil_cou)
fTPM <- readRDS(prefil_tpm)
metadata <- readRDS(metadata.rds)
metadata <- metadata %>% mutate(timetype=paste0(time,".",type)) 
# BESTIAL approach:
BEAST <- DESeqDataSetFromMatrix(countData=fmat,
                               colData=metadata,
                               design= ~ age + time + type )
keep <- rowSums(counts(BEAST)) >= 10
BEAST <- BEAST[keep,]
dsTC <- DESeq(BEAST, test="LRT", reduced= ~time+type  )
resTC <- results(dsTC)
colnames(mcols(dsTC))

d <- plotCounts(dsTC, gene=which.min(resTC$padj), intgroup="age", 
                returnData=TRUE)
d$type <- metadata[match(rownames(metadata),rownames(d)),]$type
d$time <- metadata[match(rownames(metadata),rownames(d)),]$time
pdf("gene_diff.pdf", width=8, height = 6)
ggplot(d, aes(x=age, y=count, color=age)) + 
  geom_point(position=position_jitter(w=0.1,h=0) ) + 
  scale_y_log10(breaks=c(25,100,400)) + 
  facet_grid(vars(time),vars(type)) + ggtitle("Csprs (ENSMUSG00000062783)")
dev.off()
#resLFC <- lfcShrink(dds, coef="condition_treated_vs_untreated", type="apeglm")
[1] "baseMean"            "baseVar"             "allZero"             "dispGeneEst"        
[5] "dispGeneIter"        "dispFit"             "dispersion"          "dispIter"           
[9] "dispOutlier"         "dispMAP"             "Intercept"           "time_D2_vs_D0"      
[13] "time_D4_vs_D0"       "time_D7_vs_D0"       "age_Young_vs_Old"    "SE_Intercept"       
[17] "SE_time_D2_vs_D0"    "SE_time_D4_vs_D0"    "SE_time_D7_vs_D0"    "SE_age_Young_vs_Old"
[21] "LRTStatistic"        "LRTPvalue"           "fullBetaConv"        "reducedBetaConv"    
[25] "betaIter"            "deviance"            "maxCooks"            "replace"  

