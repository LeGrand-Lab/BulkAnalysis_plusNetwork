##
# Tau specificity index calculation
# --
# johaGL 2021
#
# This script will save into different csv file:
#     housekeeping genes
#     specific genes
# in this way, separately for young and old:
  #       ids    symbol    Tau   class  whichMAX   nbMAX    day
# if several tissues exhibit max logTPM, they are all registered separated 
# by commas as single string in column whichMAX
# so 'housekeeping' and 'specific' are going to appear at the tau_class column
library(dplyr)
library(tidyverse)
library(ggplot2)
library(cowplot)
library(RColorBrewer)
library(pheatmap)
library(DESeq2)
library(ggsci) # publishing palettes
library(gridExtra)
library(reshape2)
library(ggthemes)
library(scales) # 'viridis_pal' et autres pal
library(ggrepel) # for labels to points
library("BiocParallel")
register(MulticoreParam(12)) # TODO:  set n of cores depending of available

setwd("~/BulkAnalysis_plusNetwork/")

genes_df <- read.table("data/genesinfo.csv",sep="\t",header=T)

# minimal example by hand tau specificity from Yanai:
vec <- c(0, 8, 0, 0, 0, 2, 0, 2, 0, 0, 0, 0)
num <- sum(1-(vec/max(vec)))
den <- length(vec)-1
num/den
###
# we want, for type&day,´ to obtain sample specific genes, 
# and done separately by age (old separated from young)
# import TPM matrices by day by age
# do Tau by cell type but also by celltype&day

calculateTau <- function(vec){
  if(max(vec) == 0){ # avoid zero division error
    tau.index <- 0
  }else{  
    tau.index <- (sum(1-(vec/max(vec))))/(length(vec)-1)
  }
  return(tau.index)
}

# function that parses the all days for given age, TWO different outputs:
#   1. saves the Tau dataframes to Tau/ ,.txt  separately by day
#  AND 2. returns TPM matrices into a list
saveTau.bytissue <- function(age, cts = list(), days=c("D0")){
  listofmatrices <- list() # will be returned by this function
  for (i in days){
    ktab <- read.table(paste0("data/meanTPM_",age,i,".txt"), sep='\t', header=T,
                       row.names=1)
    listofmatrices[[i]] <- ktab  # ** fill listofmatrices 
    logtab <- log(ktab+1)
    logtab[logtab <= 0] <- 0 
    tau_res <- tibble("id"=rownames(logtab))
    tau_res$symbol <- genes_df[match(rownames(logtab),genes_df$Geneid),]$symbol
    tau_res$Tau <- apply(logtab, 1, function(row) calculateTau(row))
    tau_res <- tau_res %>% mutate(class = case_when(
      Tau >= 0.8 ~ "specific",
      Tau >= 0.3 & Tau < 0.8 ~ "intermediate",
      Tau < 0.3 ~ "housekeeping"
    ))
    tissues <- colnames(logtab) # or ktab, the same
    whichMAX <- c()
    nbMAX <- c()
    for (j in 1:dim(logtab)[1]){
      roundedrow <- round(logtab[j,], 4)
      wmax <- which(roundedrow == max(roundedrow)) # position(s) in this row
      if (length(wmax)>1){
        maxtissues <- sapply(wmax, function(x) tissues[x])
        print(paste0(maxtissues, collapse=","))
        whichMAX <- c(whichMAX, paste0(maxtissues,collapse=","))
      }else{
        whichMAX <- c(whichMAX,tissues[wmax])
      }
      nbMAX <- c(nbMAX, length(wmax))
    }
    tau_res <- cbind(tau_res,whichMAX,nbMAX)
    tau_res$day <- i
    cts[[i]] <- tau_res
    write.table(bind_rows(cts), paste0("Tau/TauSpecificity_",age,i,".txt"), sep='\t',
                col.names=T, row.names=T)
  }# end for i in days
  return(listofmatrices) # ** filled by line 61
}
# here matrices go into variables : 
youngMatrices <- saveTau.bytissue("Young", days=c("D0","D2","D4","D7"))
oldMatrices <- saveTau.bytissue("Old", days =c("D0","D2","D4","D7") )

### =============== both by day and by age:
# this requires Matrices obtained above, 
# adjusts colnames to be : "M2.D7"   "FAPs.D7" "sCs.D7"  "ECs.D7" 
# instead ""M2"   "FAPs" "sCs"  "ECs" 
# and does fussion
saveTPM.bytime.tiss <- function(age,lom){  #lom : list of matrices
  dayvec <- names(lom)
  i <- 1
  fussi <- lom[[dayvec[i]]]
  newcols <- paste0(colnames(fussi),".",dayvec[i])
  colnames(fussi) <- newcols
  fussi$id <- rownames(fussi)
  for (i in 2:length(names(lom))){
    basecols <- colnames(lom[[i]])
    newcols <- paste0(basecols,".",dayvec[i])
    colnames(lom[[i]]) <- newcols
    lom[[i]]$id <- rownames(lom[[i]])
    fussi <- full_join(fussi, lom[[i]], by=c("id"))
  }
  fussi[is.na(fussi)] <- 0
  write.table(fussi, paste0("Tau/Tau_TimeType/meanTPMalldays_",age,".txt"), 
              sep='\t', col.names=T, row.names=T)
}
saveTPM.bytime.tiss("Young", youngMatrices)
saveTPM.bytime.tiss("Old", oldMatrices)

# calculate for "meanTPMalldays
tau.bytimeandtype <- function( fileprefix,age,extension){
  t <- read.table(paste0(fileprefix,age, extension), 
                  sep='\t', header=T)
  print(head(t))
  rownames(t) <- t$id
  t$id <- NULL
  log.t <- log(t+1)
  samples <- colnames(log.t)
  print(samples)
  whichMAX <- c()
  nbMAX <- c()
  for (j in 1:dim(log.t)[1]){
    roundedrow <- round(log.t[j,], 4)
    wmax <- which(roundedrow == max(roundedrow)) # position(s) in this row
    if (length(wmax)>1){
      maxtissues <- sapply(wmax, function(x) tissues[x])
      print(paste0(maxtissues, collapse=","))
      whichMAX <- c(whichMAX, paste0(maxtissues,collapse=","))
    }else{
      whichMAX <- c(whichMAX,tissues[wmax])
    }
    nbMAX <- c(nbMAX, length(wmax))
  }
  write.table(paste0("Tau/TauSpecificity_TimeType_",age,".txt"), 
              sep='\t', col.names=T, row.names=T)
}

mooo <- tau.bytimeandtype( "Tau/Tau_TimeType/meanTPMalldays_", "Young",".txt")
write.table(paste0("Tau/TauSpecificity_TimeType_",age,".txt"), 
            sep='\t', col.names=T, row.names=T)


## ====================================================================
# Filter Tau matrices : 
## ====================================================================


ages = c("Young", "Old")
days = c("D0","D2","D4","D7")
for (age in ages){
  for (day in days){
    itab <- read.table(paste0("Tau/TauSpecificity_",age,day,".txt"), sep='\t',
                       header = T)
    itab <- itab %>% filter(Tau >= 0.8) 
    if(length(itab$symbol)==unique(length(itab$symbol))){
      print("ok happy news")
    }
    print(dim(itab))
    write.table(itab, paste0("Tau/filtered/TauSpecificity_",age,day,".txt"), sep='\t',
                col.names=T, row.names=T)
  }
}





