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
library("apeglm") # BiocManager::instal("apeglm")
library("BiocParallel")
library("ashr") # install.packages("ashr")
register(MulticoreParam(12))

# using 'apeglm' for LFC shrinkage. If used in published research, please cite:
#   Zhu, A., Ibrahim, J.G., Love, M.I. (2018) Heavy-tailed prior distributions for
# sequence count data: removing the noise and preserving large differences.
# Bioinformatics. https://doi.org/10.1093/bioinformatics/bty895

# using 'ashr' for LFC shrinkage. If used in published research, please cite:
#   Stephens, M. (2016) False discovery rates: a new deal. Biostatistics, 18:2.
# https://doi.org/10.1093/biostatistics/kxw041

#
# devtools::install_github("vqv/ggbipplot")
# install.packages("glmpca")
# ** see also "signatures"

# latest : http://www.bioconductor.org/packages/release/bioc/vignettes/DESeq2/inst/doc/DESeq2.html#if-i-have-multiple-groups-should-i-run-all-together-or-split-into-pairs-of-groups
# ** cool paper: Temporal Dynamic Methods for Bulk RNA-Seq Time Series Data
# !yeah ==> https://www.biorxiv.org/content/10.1101/2020.02.01.930602v1.full
# ! yeah : https://bmcbioinformatics.biomedcentral.com/articles/10.1186/s12859-016-1403-0
# GSE : iDEA can be equally applied to bulk RNA-seq 
### >>> there exist important technical variability between replicates !

## https://bioinformatics-core-shared-training.github.io/cruk-summer-school-2018/RNASeq2018/html/06_Gene_set_testing.nb.html#fgsea
genes_df <- read.table("genesinfo.csv",sep="\t",header=T)
setwd("~/bulk_analysis/")
prefil_cou <- "prefiltered_counts.rds"
prefil_tpm <- "prefiltered_TPM.rds"
metadata.rds <- "metadata.rds"
fmat <- readRDS(prefil_cou)
fTPM <- readRDS(prefil_tpm)
metadata <- readRDS(metadata.rds)
metadata <- metadata %>% mutate(timetype=paste0(time,".",type)) 

# take away all-zeroes rows or only one count by row
fTPM <- fTPM[!rowSums(fmat)<=1,]
fmat <- fmat[!rowSums(fmat)<=1,]
# at least 3 samples with a count of 5 or higher
keep <- apply(fmat, 1, function(row) ifelse(count(row >=5)>= 3,
                                            TRUE, FALSE) )
fTPM <- fTPM[keep,]
fmat <- fmat[keep,]
# global mean and variances : all samples by old and young
df_resu.mat <- data.frame("mean.old"=rowMeans(
  fmat[,str_detect(colnames(fmat),"Old")]))
df_resu.mat$variances.old <- rowVars(
  fmat[,str_detect(colnames(fmat),"Old")]
)
df_resu.mat$mean.young = rowMeans(
  fmat[,str_detect(colnames(fmat),"Young")]
)
df_resu.mat$variances.young <- rowVars(
  fmat[,str_detect(colnames(fmat),"Young")]
)
# visualize mean count
df_resu.mat$id <- rownames(df_resu.mat)
df_resuMelt <- reshape2::melt(df_resu.mat,id=c("id"))
df_resuMelt <- df_resuMelt %>% mutate(metric= case_when(
  str_detect(variable,"mean") ~ "mean",
  str_detect(variable,"variance") ~ "variance"
) )
df_resuMelt <- df_resuMelt %>% mutate(age= case_when(
  str_detect(variable,"old") ~ "Old",
  str_detect(variable,"young") ~ "Young"
) )
ggplot(df_resuMelt, aes(age,value)) +
  geom_boxplot() + 
  geom_jitter(size=.1, alpha=.1, color="lightgray") + 
  scale_y_log10() +
  facet_wrap(~metric) + labs(title="Full dataset" ,
                             subtitle="Mean and Variance (logarithmic)")

ggplot(df_resuMelt, aes(age,value)) +
  geom_boxplot() + 
  geom_jitter(size=.1, alpha=.1, color="lightgray") + 
  scale_y_log10() +
  labs(title="Mean only")

#####
# DE analysis
# static snapshots old vs young
# by cell-type by day
#####
# filtering only the 
mergedtta <- paste0(metadata$timetype,".",metadata$age)
fixlevels <- sort(unique(mergedtta))
fixlevels <- str_replace(fixlevels, "Old", "x")
fixlevels <- str_replace(fixlevels,"Young", "Old")
fixlevels <- str_replace(fixlevels,"x", "Young")
metadata$threefx <- factor(mergedtta, 
                        levels=fixlevels,
                           ordered = F)
ddsTC <- DESeqDataSetFromMatrix(countData=fmat,
                                colData=metadata,
                                design= ~ threefx   )
(levels(colData(ddsTC)$threefx))

ddsTC <- DESeq(ddsTC, test="Wald", parallel = T )

resultsNames(ddsTC)
saveRDS(ddsTC, "DE_OldYoung_pairedDT.rds")

metagroups <- sort(unique(metadata$timetype))
# resTC_l2II <- list()
# for (i in 3:3) {  # previous run command: 
#   tmpres <- lfcShrink(ddsTC, contrast=c("threefx",
#                                         paste0(metagroups[i],".Old"),
#                                         paste0(metagroups[i],".Young")),
#                       type="ashr",
#                       parallel=T)
#   tmpres$id <- rownames(tmpres)
#   tmpres$contrast <- metagroups[i]
#   resTC_l2II[[metagroups[i]]] <- as_tibble(tmpres) %>% filter(padj<=0.7) %>%
#     arrange(padj)
# }
resTC_l2 <- list()
for (i in 1:length(metagroups)){
  tmpres <- lfcShrink(ddsTC, contrast=c("threefx",
                 paste0(metagroups[i],".Old"),
                 paste0(metagroups[i],".Young")),
                type="normal", lfcThreshold = 0.1,
                     parallel=T)
  tmpres$id <- rownames(tmpres)
  tmpres$contrast <- metagroups[i]
  resTC_l2[[metagroups[i]]] <- as_tibble(tmpres) %>% filter(padj<=0.7) %>%
    arrange(padj)
}
saveRDS(resTC_l2,"DE_OldYoung_pairedDT_BIGLIST__II.rds")


# transform this list of dataframes into single dataframe
resTC <- dplyr::bind_rows(resTC_l2)
resTC$time <- sapply(resTC$contrast, function(x) {str_split(x,'\\.')[[1]][1]})
resTC$type <- sapply(resTC$contrast, function(x) {str_split(x,'\\.')[[1]][2]})
genesymbols <- genes_df[match(resTC$id, genes_df$Geneid),]$symbol
resTC$symbols <- genesymbols
# category for genes:
resTC <- resTC %>% mutate(DEsignificant=case_when(
  padj <= 0.05 & abs(log2FoldChange) >= 1.2 ~ "FDR <= 0.05, abs(LFC)>1.2",
  padj <= 0.05 & abs(log2FoldChange) < 1.2 ~ "Not Sig",
  padj > 0.05  ~ "Not Sig"
))

write.table(resTC,"DE_table_OldYoung_pairedDT__II.csv")
###
## aesthetical modifications to get cool volcanoes
##
# n number of UP or DOWN regulated genes
DEgenUP.type <- resTC %>% 
  filter(padj <= 0.05 & lfc > 1.2) %>%
  group_by(type) %>% tally(n="UP")
DEgenDOWN.type <- resTC %>% 
  filter(padj <= 0.05 & log2FoldChange < -1.2) %>%
  group_by(type) %>% tally(n="DOWN")
DEgene.count.type <- merge(DEgenDOWN.type, DEgenUP.type, 
                                by="type")
DEgene.count.type <- paste0(DEgene.count.type$type," (",DEgene.count.type$DOWN,"|",
                            DEgene.count.type$UP,")")
names(DEgene.count.type) <- factor(c("ECs","FAPs","sCs","M1","M2","Neutro"))
#resTC$DEsignificant <- factor(resTC$significant,
#                            levels=c("FDR <= 0.0005","0.0005-0.005",
#                                     "0.005-0.05", "Not Sig"), ordered=F)
resTC$DEsignificant <- factor(resTC$DEsignificant, 
                     levels=c("FDR <= 0.05, abs(LFC)>1.2", "Not Sig"), ordered=F)
col_vir <- viridis_pal(begin=0,end=1)(10)  # nice scale, to pick from

pdf("testB.pdf",width=12,height=10)
ggplot(resTC, aes(x=log2FoldChange, y = -log10(padj))) +
  geom_point(aes(color=DEsignificant),size=.3) +
  scale_color_manual(values=c(col_vir[2],"gray")) +
  facet_grid(vars(time),vars(type),
             labeller=labeller(type=DEgene.count.type)) + 
  theme_calc() +
  theme(panel.grid.major=element_blank()) +
  geom_vline(xintercept = c(1.2,-1.2),color=col_vir[5],
             linetype="dotted", size=.5) +
  geom_text_repel(
    data= subset(resTC, padj < 0.05 & abs(log2FoldChange) > 1.2),
    aes(label=symbols, color=DEsignificant),
    size=2,
    segment.size = .1,
    force=2,
    max.overlaps=400
  ) +
  labs(title="Old vs Young across time&type",
       caption="vertical lines: ABS(log2FoldChange)=1.2 \n (DOWN | UP) regulated genes")
dev.off()


pdf("testB2.pdf",width=12,height=10)
ggplot(resTC, aes(x=log2FoldChange, y = -log10(padj))) +
  geom_point(aes(color=DEsignificant),size=.3) +
  scale_color_manual(values=c(col_vir[2],"gray")) +
  facet_grid(vars(time),vars(type),
             labeller=labeller(type=DEgene.count.type)) + 
  theme_calc() +
  theme(panel.grid.major=element_blank()) +
  geom_vline(xintercept = c(1.2,-1.2),color=col_vir[5],
             linetype="dotted", size=.5) +
  geom_text_repel(
    data= subset(resTC, padj < 0.05 & abs(log2FoldChange) > 1.2),
    aes(label=symbols, color=DEsignificant),
    size=2,
    segment.size = .1,
    force=2,
    max.overlaps=400
  ) +
  labs(title="Old vs Young across time&type",
       caption="vertical lines: ABS(log2FoldChange)=1.2 \n (DOWN | UP) regulated genes")
dev.off()

#+ theme_linedraw() +
#  theme(strip.background =element_rect(fill="white"))+
#  theme(strip.text = element_text(colour = 'black'))
# 
#
#

### recyble bin:
  # ))
  # tmpres <- results(ddsTC, contrast=c("threefx",
  #                                     paste0(metagroups[i],".Old"),
  #                                     paste0(metagroups[i],".Young")))

# resnames <- resultsNames(ddsTC)
# for (i in 1:length(metagroups)){
#   print(i)
#   for (j in 2:length(resnames)){
#     bou <- paste0("threefx_",
#                   metagroups[i],".Old_vs_", 
#                   metagroups[i],".Young")
#     if (resnames[j] == bou){
#       print("ok")
#       print(bou)
#     }
#   }
# }

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

##### another test

ddsTC <- DESeqDataSetFromMatrix(countData = fmat,
                                colData = metadata,
                                design = ~ type + age + time + age:time  )
ddsTC <- DESeq(ddsTC, test="LRT", reduced = ~ type + time)
resTC <- results(ddsTC)
ggthing <- plotCounts(ddsTC, which.min(resTC$padj), 
                   intgroup = c("time","age"), returnData = TRUE)
ggthing$type <- metadata[match(rownames(metadata),rownames(ggthing)),]$type
ggplot(ggthing,
       aes(x = time, y = count, color = age, group = age)) + 
  geom_point() + stat_summary(fun=mean, geom="line") +
  scale_y_log10() + facet_grid(~type)

resTC$index <- seq(1:dim(resTC)[1])
ggthing2 <- plotCounts(ddsTC, resTC['ENSMUSG00000019932',]$index, 
                      intgroup = c("time","age"), returnData = TRUE)
ggthing2$type <- metadata[match(rownames(metadata),rownames(ggthing2)),]$type
ggplot(ggthing2,
       aes(x = time, y = count, color = age, group = age)) + 
  geom_point() + stat_summary(fun=mean, geom="line") +
  scale_y_log10() + facet_grid(~type)


savedtc <- resTC
View(as.data.frame(savedtc[order(savedtc$padj),]))
###

ddsTC <- DESeqDataSetFromMatrix(countData = fmat,
                                colData = metadata,
                                design = ~ type + time + age  )
ddsTC <- DESeq(ddsTC, test="LRT", reduced= ~type + time)
resTC <- results(ddsTC)
resultsNames(ddsTC)


ggthing3 <- plotCounts(ddsTC, which.min(resTC$padj), 
                      intgroup = c("time","age", "type"), returnData = TRUE)
ggplot(ggthing3,
       aes(x = time, y = count, color = age, group = age)) + 
  geom_point() + stat_summary(fun=mean, geom="line") +
  scale_y_log10() + facet_grid(~type)
