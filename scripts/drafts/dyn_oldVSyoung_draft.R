## 
# Detect Differential Expression
# --
# Joha GL
##
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
library(reshape2)
library(ggthemes)
library(scales) # 'viridis_pal' et autres pal
library(ggrepel) # for labels to points
library("apeglm") # BiocManager::install("apeglm")
library("BiocParallel")
library("ashr") # install.packages("ashr")
register(MulticoreParam(4)) # TODO:  set n of cores depending of available
library(pheatmap)
setwd("~/bulk_analysis/")
resdir = "fromhome/"
genes_df <- read.table("genesinfo.csv",sep="\t",header=T)

prefil_cou <- "prefiltered_counts.rds"
prefil_tpm <- "prefiltered_TPM.rds"
metadata.rds <- "metadata.rds"
fmat <- readRDS(prefil_cou)
fTPM <- readRDS(prefil_tpm)
metadata <- readRDS(metadata.rds)
metadata <- metadata %>% mutate(timetype=paste0(time,".",type)) 

keep <- apply(fmat, 1, function(row) ifelse(count(row >=5)>= 3, TRUE, FALSE) )
fTPM <- fTPM[keep,]
fmat <- fmat[keep,]

# divide into celltypes, 'dynamic' exp here only M2
ct <- "M2"
mat.ct <- fmat[,str_detect(colnames(fmat), ct)]
meta.ct <- metadata %>% filter(str_detect(rownames(metadata),ct))

x.keep <- apply(mat.ct, 1, function(row) ifelse(sum(row >=5)>= 3, TRUE, FALSE))
mat.ct <- mat.ct[x.keep,]
ds.o <- DESeqDataSetFromMatrix(countData = mat.ct,
                               colData = meta.ct,
                               design = ~ age + time + age:time)
colData(ds.o)$age <- relevel(colData(ds.o)$age, ref="Young")
ddsTC <- DESeq(ds.o, test="LRT", reduced = ~ age + time)
resTC <- results(ddsTC, parallel=T)
resTC$symbol <- mcols(ddsTC)$symbol
head(resTC[order(resTC$padj),], 4)

ggco <- plotCounts(ddsTC, which.min(resTC$padj),
                   intgroup = c("time", "age"),
                   returnData = T)
#ggco$time <- as.character(ggco$time)
ggplot(ggco,
       aes(x=time, y=count, color= age,
           group=age)) +
  geom_point() + stat_summary(fun.y=mean, geom="line") +
  scale_y_log10()

tdf <- as_tibble(resTC) %>% mutate(id=rownames(resTC)) %>%
  filter(padj < 0.05 & abs(log2FoldChange) >= 1.2)
# there are 7 genes significatively dyn exp but to get
# information about which one goes up/down at D2/D4/D7???
resultsNames(ddsTC)
betas <- coef(ddsTC)

topGenes <- head(order(resTC$padj),20)
mat <- betas[topGenes, -c(1,2)]
thr <- 3 
mat[mat < -thr] <- -thr
mat[mat > thr] <- thr
pheatmap(mat, breaks=seq(from=-thr, to=thr, length=101),
         cluster_col=FALSE)
## ces contrastes:  ageOldtimeD7, ageOldtimeD4 timeD7vsD2 timeD4vsD2 sont confus !!!
# TODO : faire des meilleures contrastes... maybe
# TODO: mettre les gene symbols