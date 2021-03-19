## 
# Detect Markers at D7 in young
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
library(gridExtra)
library(reshape2)
library(ggthemes)
library(scales) # 'viridis_pal' et autres pal
library(ggrepel) # for labels to points
library("apeglm") # BiocManager::install("apeglm")
library("BiocParallel")
library("ashr") # install.packages("ashr")
register(MulticoreParam(12)) # TODO:  set n of cores depending of available

setwd("~/bulk_analysis/")
odir = "signaturestypes/"
genes_df <- read.table("genesinfo.csv",sep="\t",header=T)

symbols <- 'data/genesinfo.csv'
symb_df <- read.table(symbols, header=T)
prefil_cou <- "data/prefiltered_counts.rds"
prefil_tpm <- "data/prefiltered_TPM.rds"
metadata.rds <- "data/metadata.rds"
fmat <- readRDS(prefil_cou)
fTPM <- readRDS(prefil_tpm)
metadata <- readRDS(metadata.rds)
metadata <- metadata %>% mutate(timetype=paste0(time,".",type)) 

keep <- apply(fmat, 1, function(row) ifelse(count(row >=5)>= 3, TRUE, FALSE) )
fTPM <- fTPM[keep,]
fmat <- fmat[keep,]

# get young only
xmat <- fmat[,str_detect(colnames(fmat),"Young")]
xmeta <- metadata[str_detect(rownames(metadata),"Young"),]
rm(metadata, fmat)

# desired contrasts: 
#      -"baseline" for faps ec and sat is D0
#      -no baseline for M1 M2 , take D7 as "homeostasis"
#      -Neutrophyls are only in D2
ALLpairs_day <- function(day, mat, metadata){
  sortie <- data.frame("baseMean"= double(), "log2FoldChange"= double(), 
                       "lfcSE"= double(),  "stat" = double(), 
                       "pvalue" = double(), "padj" = double(),
          "id" = character(),  "query" = character(), "ref" = character(),
                    "day" = character())
  xmat <- mat[,str_detect(colnames(mat),day)]
  xmeta <- metadata[str_detect(rownames(metadata),day),]
  
  x.keep <- apply(xmat, 1, function(row) ifelse(sum(row >=5)>= 3, TRUE, FALSE))
  xmat <- xmat[x.keep,]
  ds.o <- DESeqDataSetFromMatrix(countData = xmat,
                                 colData = xmeta,
                                 design = ~ type)
  # initially set satellite as reference (because 'stem' type)
  nostem <- setdiff( c(levels(colData(ds.o)$type)), c("sCs"))
  colData(ds.o)$type <- factor(colData(ds.o)$type, 
                               levels = c(c("sCs"), c(nostem) ) )
  dds <- DESeq(ds.o, parallel=T)
  alltyp <- as.character(levels(colData(ds.o)$type))
  for(i in alltyp){
    therest = setdiff(alltyp,i)
    print(therest)
    for (j in therest){
      rds <- results(dds, contrast = c("type", j, i) ,
                     parallel=T)
      tmpdf <- as_tibble(rds) %>% mutate(id=rownames(rds)) %>%
        filter(padj < 0.05 & abs(log2FoldChange) >= 1.2)
      tmpdf$query <- j 
      tmpdf$ref <- i
      tmpdf$day <- day
      sortie <- rbind(sortie,tmpdf)
      print(dim(sortie))
    }
  }
  return(sortie)
}
# NOTE: Not possible to use (saddly):  list to contrast and 'listValues'              
# i.e. contrasting the first group vs the average of the other groups
# that's why I made paired contrasts, then selecting most important markers
# (complex processing of output dataframe)
day <- "D7"
D7.paired.Y = ALLpairs_day(day, xmat, xmeta)

write.table(D7.paired.Y, file=paste0(odir,"D7_allcontrasts_Y.csv"), col.names = T,
            row.names = F)

#### open saved table
D7.paired.Y <- read.table(paste0(odir,"D7_allcontrasts_Y.csv"),header=T)
dim(D7.paired.Y)

## make groups of genes by query celltype:
D7.paired.Y$query <- factor(D7.paired.Y$query )
(D7.paired.Y$query) 
D7.paired.Y$symbol <- symb_df[match(D7.paired.Y$id, symb_df$Geneid,),]$symbol

test_dedup <-  D7.paired.Y %>% group_by(query) %>%
  dplyr::distinct(symbol, query, .keep_all=T)  # just keeps first occurrence
#filter this enormous dataframe
test_dedup <- test_dedup %>% filter(log2FoldChange > 1.2 & padj < 1e-20) %>%
  group_by(query) %>% 
  arrange(padj, by_group=T) 
test_dedup <- test_dedup %>% group_by(query) %>% arrange(log2FoldChange) %>%
  slice_max(log2FoldChange, n=700)


## extract a minimal subset to do the interactome:
minitest <- test %>% group_by(query) %>%  slice_max(log2FoldChange, n=20)
## i found italk, which allows to do differential exp between old-young
#  so il will extract matrix by datapoint containing both conditions
