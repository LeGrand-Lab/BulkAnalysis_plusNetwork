## 
# Detect Differential Expression at each time-point, by cellType
# --
# Joha GL
##
library(dplyr)
library(tidyverse)
library(ggplot2)
library(cowplot)
library(forcats)
library(RColorBrewer)
library(pheatmap)
library(DESeq2)
library(ggsci) # publishing palettes
library(gridExtra)
library(reshape2)
library(ggthemes)
library(scales) # 'viridis_pal' et autres pal
library(ggrepel) # for labels to points
library("apeglm") # BiocManager::install("apeglm")
library("BiocParallel")
library("ashr") # install.packages("ashr")
register(MulticoreParam(4)) # TODO:  set n of cores depending of available
# devtools::install_github("vqv/ggbipplot")
# install.packages("glmpca")

# using 'apeglm' for LFC shrinkage. If used in published research, please cite:
#   Zhu, A., Ibrahim, J.G., Love, M.I. (2018) Heavy-tailed prior distributions for
# sequence count data: removing the noise and preserving large differences.
# Bioinformatics. https://doi.org/10.1093/bioinformatics/bty895

# using 'ashr' for LFC shrinkage. If used in published research, please cite:
#   Stephens, M. (2016) False discovery rates: a new deal. Biostatistics, 18:2.
# https://doi.org/10.1093/biostatistics/kxw041

# latest : http://www.bioconductor.org/packages/release/bioc/vignettes/DESeq2/inst/doc/DESeq2.html#if-i-have-multiple-groups-should-i-run-all-together-or-split-into-pairs-of-groups
# ** cool paper: Temporal Dynamic Methods for Bulk RNA-Seq Time Series Data
# !yeah ==> https://www.biorxiv.org/content/10.1101/2020.02.01.930602v1.full
# ! yeah : https://bmcbioinformatics.biomedcentral.com/articles/10.1186/s12859-016-1403-0
# GSE : iDEA can be equally applied to bulk RNA-seq 
### >>> there exist important technical variability between replicates !

## https://bioinformatics-core-shared-training.github.io/cruk-summer-school-2018/RNASeq2018/html/06_Gene_set_testing.nb.html#fgsea
setwd("~/BulkAnalysis_plusNetwork/")
odir = "exam_INTER_conditions/static/"
shot_t = "shot_dataframe_unfiltered.csv"
shot_t_f = "shot_dataframe.csv"

genes_df <- read.table("data/genesinfo.csv", sep="\t", header=T)

prefil_cou <- "data/prefiltered_counts.rds"
metadata.rds <- "data/metadata.rds"
fmat <- readRDS(prefil_cou)
metadata <- readRDS(metadata.rds)
metadata <- metadata %>% mutate(timetype=paste0(time,".",type)) 
keep <- apply(fmat, 1, function(row) ifelse(count(row >=5)>= 3, TRUE, FALSE) )
fmat <- fmat[keep,]

#####
# DE analysis
# static snapshots old vs young
# by cell-type by day
#####

##  with 'subset' s   on same original object:
#  ===========================================================================
ds.o <- DESeqDataSetFromMatrix(countData = fmat,
                               colData = metadata,
                               design= ~ age)
# set "Young" as reference
colData(ds.o)$age <- factor(colData(ds.o)$age,
                            levels=c("Young", "Old"))

time.type <- sort(unique(colData(ds.o)$timetype))

sortie <- data.frame("baseMean"= double(), "log2FoldChange"= double(), 
                     "lfcSE"= double(),  "stat" = double(), 
                     "pvalue" = double(), "padj" = double(),
                     "id" = character(), "day" = character(), "type"=character())

for (k in time.type){
  time <- str_split(k, "\\.")[[1]][1]
  type <- str_split(k, "\\.")[[1]][2]
  print(paste(time,type))
  # do subsets as needed:
  ds.sub.o <- ds.o[,ds.o$type == type & ds.o$time == time]
  ds.sub <- DESeq(ds.sub.o, test="Wald", full=~age)
  deres <- lfcShrink(ds.sub, coef="age_Old_vs_Young", type="apeglm",
                     parallel=T)
  # run thing:
  deres$id <- rownames(deres)
  deres <- as_tibble(deres) %>% filter(padj <= 0.5 & abs(log2FoldChange) >= 0.1)
  deres$day = time
  deres$type = type
  sortie <- rbind(sortie,deres)
}
# save these results
write.table(sortie, paste0(odir, shot_t), sep="\t", col.names=T, row.names = F)

# save significant extract (~1000 genes):
write.table(sortie %>% filter(padj<0.05 & abs(log2FoldChange) >= 1.2),
            paste0(odir,shot_t_f), sep="\t", col.names=T, row.names = F)
# complete table with gene_symbols:
sortie <- read.table(paste0(odir,shot_t_f), sep="\t", header=T)

rownames(genes_df) <- genes_df$Geneid
symvec <- c()
for (i in sortie$id){
  symvec <- c(symvec, genes_df[i,]$symbol)
}
sortie$symbol <- symvec
write.table(sortie %>% filter(padj<0.05 & abs(log2FoldChange) >= 1.2),
            paste0(odir,shot_t_f), sep="\t", col.names=T, row.names = F)
####  downstream
# how many UP?
upsortie <- sortie %>% filter(log2FoldChange > 0.05)
(dim(upsortie))
View(upsortie)

# DOWN?
dosortie <- sortie %>% filter(log2FoldChange < 0)

# ** VOLCANOES ** :
sortie <- read.table(paste0(odir, shot_t),sep="\t", header=T)
# set aesthetics data:
sortie <- sortie %>% mutate(DEsignificant=case_when(
  padj < 1e-10 & abs(log2FoldChange) >= 1.2 ~ "FDR < 1e-10" ,
  padj <= 0.05 & padj > 1e-10 & abs(log2FoldChange) >= 1.2 ~ "FDR <= 0.05 ",
  padj > 0.05 | abs(log2FoldChange) < 1.2 ~ "Not Sig"
))
sortie$symbol <- genes_df[match(sortie$id, genes_df$Geneid),]$symbol
# n number of UP or DOWN regulated genes
DEgenUP.type <- sortie %>% 
  filter(padj <= 0.05 & log2FoldChange > 1.2) %>%
  group_by(type) %>% tally(n="UP")
DEgenDOWN.type <- sortie %>% 
  filter(padj <= 0.05 & log2FoldChange < -1.2) %>%
  group_by(type) %>% tally(n="DOWN")
DEgene.count.type <- merge(DEgenDOWN.type, DEgenUP.type, by="type")
nbDE.vec <- paste0(DEgene.count.type$type," (",DEgene.count.type$DOWN,"|",
                            DEgene.count.type$UP,")")
names(nbDE.vec) <- factor(DEgene.count.type$type)


col_vir <- viridis_pal(begin=0,end=1)(10)  # nice scale, to pick from

# print stuff:
pdf(paste0(odir,"volcano_snapshot.pdf"), width=14, height=10)
ggplot(sortie, aes(x=log2FoldChange, y = -log10(padj+1e-20))) +
  geom_point(aes(color=DEsignificant),size=.3) +
  scale_color_manual(values=c(col_vir[7], col_vir[3],"gray")) +
  facet_grid(vars(day),vars(type),
             labeller=labeller(type=nbDE.vec)) + 
  theme_calc() +
  theme(panel.grid.major=element_blank()) +
  geom_vline(xintercept = c(1.2,-1.2),color=col_vir[5],
             linetype="dashed", size=.2) +
  geom_text_repel(
    data= subset(sortie, padj < 0.0005 & 
                   abs(log2FoldChange) > 1.2),
    aes(label=symbol, color=DEsignificant),
    size=2,
    segment.size = .1,
    force=2, force_pull = 2,
    max.overlaps=15
  ) +
  labs(title="Old vs Young across time&type",
       caption="vertical lines: ABS(log2FoldChange)=1.2 \n 
       (DOWN | UP) regulated genes \n
       labels only for genes padj < 0.0005")
dev.off()


