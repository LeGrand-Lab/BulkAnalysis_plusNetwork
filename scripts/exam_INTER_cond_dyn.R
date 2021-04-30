## 
# Dynamics of diff expression (in time)
# Each cell type separately (Neutro EXCLUDED as single-time-point measured)
# --
# Joha GL
##
library(dplyr)
library(tidyverse)
library(RColorBrewer)
library(DESeq2)
library(cowplot)
library(reshape2)
library("BiocParallel")
register(MulticoreParam(4)) # TODO:  set n of cores depending of available


setwd("~/BulkAnalysis_plusNetwork/")
prefil_cou <- "data/prefiltered_counts.rds"
metadata.rds <- "data/metadata.rds"
genes_df <- read.table("data/genesinfo.csv",sep="\t",header=T)

resdir = "exam_INTER_conditions/dynamic/"
fmat <- readRDS(prefil_cou)
metadata <- readRDS(metadata.rds)
metadata <- metadata %>% mutate(typetimeage=paste0(type,".",time,".",age)) 

# verify batches correspond to cell type:
batchesinfo <- read.table("data/batchesinfo.csv",sep="\t",header=T)
head(batchesinfo)
metadata$batch = batchesinfo[match(metadata$sample, batchesinfo$sample),]$batch
type_batch = metadata %>% select(type,batch) %>% distinct(type,batch, .keep_all=T)
type_batch
# Old.M1.D2_1         M1    B1
# Old.M2.D7_3         M2    B1
# Young.FAPs.D2_1   FAPs    B2
# Young.sCs.D4_2     sCs    B3
# Young.ECs.D4_2     ECs    B4
# Old.Neutro.D2_2 Neutro    B4
##  set which rows to keep :
keep <- apply(fmat, 1, function(row) ifelse(count(row >=5)>= 3, TRUE, FALSE) )
fmat <- fmat[keep,]

plotDispersionsTissues <- function( fmat, metadata, tisscol = "type",
                                    outfilename, nbco, nbro){
  tissues = sort(unique(metadata[[tisscol]]))
  pdf(outfilename, paper="a4")
  par(mfrow=c(nbro,nbco))
  for (ct in tissues){
    mat.ct <- fmat[,str_detect(colnames(fmat), ct)]
    meta.ct <- metadata %>% filter(str_detect(rownames(metadata),ct))
    x.keep <- apply(mat.ct, 1, function(row) ifelse(sum(row >=5)>= 3, TRUE, FALSE))
    mat.ct <- mat.ct[x.keep,]
    ds.o <- DESeqDataSetFromMatrix(mat.ct, meta.ct, design = ~ age)
    ds.o$age <- relevel(ds.o$age, ref="Young")
    d <- DESeq2::estimateSizeFactors(ds.o)
    d <- DESeq2::estimateDispersions(d)
    plotDispEsts(d, cex=0.6, genecol=rgb(0.1,0.1,0.1,0.3), fitcol="gold", 
                 finalcol=rgb(0.1,0.8,0.7,0.6), main=ct)
  } # end for
  dev.off()
  par(mfrow=c(1,1))
}
# plotDispersionsTissues(fmat, metadata, "type", 
#                        outfilename=paste0(resdir,"dispersionsPlot.pdf"), 2, 3)

# use runInteracCustom when equal or more than 3 time-points available
# this function recalculates FDR (repadj) but filters  on padj
# NOTE : padjcutoff is loose by default
runInteracCustom <- function(ct, fmat, metadata, resdir, 
                             lfccutoff=1.2, padjcutoff=1){
  mat.ct <- fmat[,str_detect(colnames(fmat), ct)]
  meta.ct <- metadata %>% filter(str_detect(rownames(metadata),ct))
  x.keep <- apply(mat.ct, 1, function(row) ifelse(sum(row >=5)>= 3, TRUE, FALSE))
  mat.ct <- mat.ct[x.keep,]
  ds.o <- DESeqDataSetFromMatrix(mat.ct, meta.ct, design = ~ age + time + age:time)
  ds.o$age <- relevel(ds.o$age, ref="Young")
  tipoints <- as.character(sort(unique(ds.o$time)))
  conseclist <- list()
  for (i in 1:(length(tipoints)-1)){
    consec = c(tipoints[i+1],tipoints[i])
    print(consec)
    ds.consec <- ds.o[,as.character(ds.o$time) %in% consec] # subset object
    #update levels
    ds.consec$time <- factor(ds.consec$time, levels=rev(consec)) 
    ds.consec$age <- relevel(ds.consec$age, ref="Young")
    d.r <- DESeq(ds.consec, full=~age + time + age:time)
    print("ended DE step, picking interaction contrast")
    # the interaction term, answering: is the time effect *different* across ages ?
    mycontrast <- paste0("ageOld.time", consec[1])
    tint = results(d.r, name = mycontrast )
    tint$id <- rownames(tint)
    tint <- as_tibble(tint) %>% mutate(contrast=paste0(consec,collapse="_vs_"),
                                       celltype = ct ) 
    conseclist[[paste0(consec,collapse="_vs_")]] <- tint
  }
  tmpprep = bind_rows(conseclist)
  pvals = tmpprep$pvalue
  tmpprep$repadj = p.adjust(pvals, method = "BH", n = length(pvals))
  write.table(tmpprep %>% filter(abs(log2FoldChange) >= lfccutoff &
                                padj <= padjcutoff), # see non significant also
              file=paste0(resdir, ct, "_INTERagetime.csv"),
              sep='\t', col.names=T, row.names = F)
  print(paste("saved csv file for INTERaction, ",ct))
  return("done")
}

# use this _simple function if only 2 time-points available : 
runInteracCustom_simple <- function(ct, fmat, metadata, resdir, 
                             lfccutoff=1.2, padjcutoff=1){
  mat.ct <- fmat[,str_detect(colnames(fmat), ct)]
  meta.ct <- metadata %>% filter(str_detect(rownames(metadata),ct))
  x.keep <- apply(mat.ct, 1, function(row) ifelse(sum(row >=5)>= 3, TRUE, FALSE))
  mat.ct <- mat.ct[x.keep,]
  ds.o <- DESeqDataSetFromMatrix(mat.ct, meta.ct, design = ~ age + time + age:time)
  ds.o$age <- relevel(ds.o$age, ref="Young")
  timeordered <- sort(unique(ds.o$time))
  ds.o$time <- factor(ds.o$time, levels=timeordered)
  d.r <- DESeq(ds.o, full=~age + time + age:time)
  print("ended DE step, picking interaction contrast")
  # the interaction term, answering: is the time effect *different* across ages ?
  mycontrast <- paste0("ageOld.time", timeordered[2])
  tint <- results(d.r, name=mycontrast)
  print(head(tint))
  tint$id <- rownames(tint)
  tint <- as_tibble(tint) %>% 
    mutate(contrast=paste0(rev(timeordered),collapse="_vs_"),
                                     celltype = ct ) 
  write.table(tint %>% filter(abs(log2FoldChange) >= lfccutoff &
                     padj <= padjcutoff), # see non significant also
              file=paste0(resdir, ct, "_INTERagetime.csv"),
              sep='\t', col.names=T, row.names = F)
  print(paste("saved csv file for INTERaction, ",ct))
  return("done")
}


# use only Tissues having more than one timepoint:
# let loose cutoff for padj (as default ) because will undergo GSEA:
runInteracCustom("FAPs", fmat, metadata, resdir=resdir)
runInteracCustom("ECs", fmat, metadata, resdir=resdir)
runInteracCustom_simple("M1", fmat, metadata, resdir=resdir, lfccutoff = 0.1)
runInteracCustom("M2", fmat, metadata, resdir=resdir)
runInteracCustom("sCs", fmat, metadata, resdir=resdir)


## end  execution


# APPENDIX
# ===========================================================================
# set this to TRUE if an illustrative exampel desired :
neededIllustrTest = FALSE
if (neededIllustrTest){
  print("this example requires fmat and metadata objects loaded as above")
  ct = "FAPs"# *******
  mat.ct <- fmat[,str_detect(colnames(fmat), ct)]
  meta.ct <- metadata %>% filter(str_detect(rownames(metadata),ct))
  x.keep <- apply(mat.ct, 1, function(row) ifelse(sum(row >=5)>= 3, TRUE, FALSE))
  mat.ct <- mat.ct[x.keep,]
  ds.o <- DESeqDataSetFromMatrix(mat.ct, meta.ct, design = ~ age + time + age:time)
  ds.o$age <- relevel(ds.o$age, ref="Young")
  d.demo <- DESeq(ds.o)
  # time efect for reference age (main effect)
  refd = results(d.demo, contrast = c("time","D2","D0"))
  # time effect for age Old: the extra time effect in Old compared to ref.
  comd = results(d.demo, list(c("time_D2_vs_D0", "ageOld.timeD2")))
  # by hand calculate deltas
  if(all(rownames(refd)==rownames(comd))){
    demodf = data.frame(id=rownames(refd),
                        deltas =  comd$log2FoldChange - refd$log2FoldChange)}
  # Compare my deltas with true DESeq2 interaction:
  # the interaction term do test: is the time effect *different* across ages ?
  # this means :  [(OldD2vsOldD0)vs(YoungD2vsYoungD0)]
  tint = results(d.demo, name = "ageOld.timeD2")  
  tint$id <- rownames(tint)
  demodf <- inner_join(as_tibble(tint),demodf,by=c("id"))
  View(demodf %>% filter(padj <= 0.05 & abs(log2FoldChange)>1.5))
  # fantastic ! my calculated deltas are VERY close to interaction lfc found ! 
  print(max(demodf$log2FoldChange - demodf$deltas))
}
# ===========================================================================
#library(reticulate)
#use_condaenv(condaenv="~/prog_bio/anaconda3/envs/MYVENV",
#             conda="~/prog_bio/anaconda3/bin/conda")
