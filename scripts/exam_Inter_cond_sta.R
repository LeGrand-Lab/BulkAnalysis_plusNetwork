## 
# Detect Differential Expression at each time-point, by cellType
# output: all results into 'exam_INTER_conditions/static/'
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

setwd("~/BulkAnalysis_plusNetwork/")
odir = "exam_INTER_conditions/static/"

shot_sf = "shot_dataframe_softfilter.csv"  # soft filters are : padj 0.5 and lfc 0.1 
shot_fi = "shot_dataframe_finalfilter.csv" # final filter applied padj0.05 lfc>= 1.2

genes_df <- read.table("data/genesinfo.csv", sep="\t", header=T)
volcanoneeded <- FALSE
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

sortie_full <- data.frame("baseMean"= double(), "log2FoldChange"= double(), 
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
  deres$id <- rownames(deres)
  deres$day = time
  deres$type = type
  # pile full output
  sortie_full <- rbind(sortie_full, as_tibble(deres))
  deres_f <- as_tibble(deres) %>% filter(padj <= 0.5 & abs(log2FoldChange) >= 0.1)
  out_softfilt <- rbind(out_softfilt,deres_f)
}
# complete tables with gene_symbols:
sortie_full$symbol <- genes_df[match(sortie_full$id, genes_df$Geneid),]$symbol
# save these results
saveRDS(sortie_full, paste0(odir, "shot_rds_full.rds")) 
#write.table(sortie_full, paste0(odir, "shot_dataframe_full.csv"), sep="\t", col.names=T, row.names = F) 
out_softfilt <- sortie_full %>% filter(padj <= 0.5 & abs(log2FoldChange) >= 0.1)

write.table(out_softfilt, paste0(odir, shot_sf), 
            sep="\t", col.names=T, row.names = F) 
# save significant extract (~1000 genes):
significant_df = out_softfilt %>% filter(padj <= 0.05 & abs(log2FoldChange) >= 1.2)
write.table(significant_df, paste0(odir,shot_fi), sep="\t", col.names=T, row.names = F)
####  downstream
# how many UP
upsortie <- out_softfilt %>% filter(log2FoldChange > 0.05)
(dim(upsortie))
View(upsortie)
# how many DOWN
dosortie <- out_softfilt %>% filter(log2FoldChange < 0)

# ** VOLCANOES ** :
if (volcanoneeded){
  out_softfilt <- read.table(paste0(odir, shot_sf),sep="\t", header=T)
  # set aesthetics data:
  out_softfilt <- out_softfilt %>% mutate(DEsignificant=case_when(
    padj < 1e-10 & abs(log2FoldChange) >= 1.2 ~ "FDR < 1e-10" ,
    padj <= 0.05 & padj > 1e-10 & abs(log2FoldChange) >= 1.2 ~ "FDR <= 0.05 ",
    padj > 0.05 | abs(log2FoldChange) < 1.2 ~ "Not Sig"
  ))
  out_softfilt$symbol <- genes_df[match(out_softfilt$id, genes_df$Geneid),]$symbol
  # n number of UP or DOWN regulated genes
  DEgenUP.type <- out_softfilt %>% 
    filter(padj <= 0.05 & log2FoldChange > 1.2) %>%
    group_by(type) %>% tally(n="UP")
  DEgenDOWN.type <- out_softfilt %>% 
    filter(padj <= 0.05 & log2FoldChange < -1.2) %>%
    group_by(type) %>% tally(n="DOWN")
  DEgene.count.type <- merge(DEgenDOWN.type, DEgenUP.type, by="type")
  nbDE.vec <- paste0(DEgene.count.type$type," (",DEgene.count.type$DOWN,"|",
                     DEgene.count.type$UP,")")
  names(nbDE.vec) <- factor(DEgene.count.type$type)
  
    col_vir <- viridis_pal(begin=0,end=1)(10)  # nice scale, to pick from
  
  # print stuff:
  pdf(paste0(odir,"volcano_snapshot.pdf"), width=14, height=10)
  ggplot(out_softfilt, aes(x=log2FoldChange, y = -log10(padj+1e-20))) +
    geom_point(aes(color=DEsignificant),size=.3) +
    scale_color_manual(values=c(col_vir[7], col_vir[3],"gray")) +
    facet_grid(vars(day),vars(type),
               labeller=labeller(type=nbDE.vec)) + 
    theme_calc() +
    theme(panel.grid.major=element_blank()) +
    geom_vline(xintercept = c(1.2,-1.2),color=col_vir[5],
               linetype="dashed", size=.2) +
    geom_text_repel(
      data= subset(out_softfilt, padj < 0.0005 & 
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
}


# ============================================================================
# END
#
# final notes:
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

#rownames(genes_df) <- genes_df$Geneid
#symvec <- c()
#for (i in sortie_full$id){
#  symvec <- c(symvec, genes_df[i,]$symbol)
#}
#sortie_full$symbol <- symvec
