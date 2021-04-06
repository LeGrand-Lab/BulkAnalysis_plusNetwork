##
# Tau specificity index calculation
# --
# johaGL 2021
#
# This script will save into different txt and xlsx files:
#     housekeeping genes
#     specific genes
# in this way, separately for young and old:
#       ids    symbol    Tau   class  whichMAX   nbMAX    day
# if several tissues exhibit max logTPM, they are all registered separated 
# by commas as single string in column whichMAX
# so 'housekeeping' and 'specific' are going to appear at the tau_class column
library(dplyr)
library(tidyverse)
library(openxlsx)
library(ggplot2)
library(cowplot)
library(RColorBrewer)
library(ggsci) # publishing palettes
library(gridExtra)
library(reshape2)
library(ggthemes)
library(scales) # 'viridis_pal' et autres pal
library(ggrepel) # for labels to points
library("BiocParallel")
register(MulticoreParam(8)) # TODO:  set n of cores depending of available

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

# function that parses the all days for given age:
#   output: saves the Tau dataframes to Tau/ ,.txt  separately by day
saveTau.bytissue <- function(age,  days){
  for (i in days){
    ktab <- read.table(paste0("data/meanTPM_",age,i,".txt"), sep='\t', header=T,
                       row.names=1) 
    logtab <- log10(ktab+1) 
    q.cutoff <- quantile(unlist(logtab),0.10)
    keep <- apply(logtab, 1, function(x) sum(x >= 0) == length(x) &
                    sum(x > q.cutoff) >= 1)  #and at least one over this value
    logtab <- logtab[keep,] 
    print(dim(logtab))
    tau_res <- tibble("id"=rownames(logtab))
    tau_res$symbol <- genes_df[match(rownames(logtab),genes_df$Geneid),]$symbol
    tau_res$Tau <- apply(logtab, 1, function(row) calculateTau(row))
    tau_res <- tau_res %>% mutate(class = case_when(
      Tau >= 0.8 ~ "specific",
      Tau >= 0.5 & Tau < 0.8 ~ "intermediate",
      Tau < 0.5 ~ "housekeeping"
    ))
    tissues <- colnames(logtab) # or ktab, the same
    whichMAX <- c()
    nbMAX <- c()
    maxlog10TPM <- c()
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
      maxlog10TPM <- c(maxlog10TPM, max(logtab[j,]))
    }
    tau_res <- cbind(tau_res,whichMAX,nbMAX, maxlog10TPM)
    tau_res$day <- i
    write.table(tau_res, paste0("Tau/TauSpecificity_",age,i,".txt"), sep='\t',
                col.names=T, row.names=T)
  }# end for i in days
  return("tau calculated & txt files saved into 'Tau/'")
}
#print( saveTau.bytissue("Young", days=c("D0","D2","D4","D7")) )
#print( saveTau.bytissue("Old",  days=c("D0","D2","D4","D7")) ) 

## ====================================================================
# Filter generated Tau matrices : 
## ====================================================================
ages = c("Young", "Old")
days = c("D0","D2","D4","D7")
for (age in ages){
  for (day in days){
    itab <- read.table(paste0("Tau/TauSpecificity_",age,day,".txt"), sep='\t',
                       header = T)
    hist(itab$Tau, main=paste0(age,":", day))
    itab <- itab %>% filter(class == "specific") 
    if(length(itab$symbol)==unique(length(itab$symbol))){
      print("ok happy news")
    }
    print(dim(itab))
    #write.table(itab, paste0("Tau/filtered/TauSpecificity_",age,day,".txt"), sep='\t',
    #             col.names=T, row.names=T)
    # prepare xls file
    xfile = paste0("Tau/filtered/TauSpecificity_",age,day,".xlsx")
    wb <- createWorkbook()
    typesspe = unique(itab$whichMAX)
    print(typesspe)
    for (ty in typesspe){
      addWorksheet(wb, ty )
      writeDataTable(wb, sheet = ty , 
                     x = itab %>% filter(whichMAX==ty & 
                                           maxlog10TPM > 0.5),
                     colNames = TRUE, rowNames = F)
    }
    saveWorkbook(wb, xfile, overwrite = TRUE)
    #openXL("...xlsx")
  }
}
## ====================================================================

# ########################## EXPERIMENTAL
### =============== both by day and by age:
# this requires Matrices obtained above, 
# adjusts colnames to be : "M2.D7"   "FAPs.D7" "sCs.D7"  "ECs.D7" 
# instead ""M2"   "FAPs" "sCs"  "ECs" 
# and does fussion to obtain 2 dataframes: one for age
saveTPM.bytime.tiss <- function(age,cts = list(), days=c("D0")){ 
  lom <- list() #lom : list of matrices
  for (i in days){
    ktab <- read.table(paste0("data/meanTPM_",age,i,".txt"), sep='\t', header=T,
                       row.names=1) 
    lom[[i]] <- ktab  # ** fill listofmatrices 
  }
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

saveTPM.bytime.tiss("Young", days=c("D0","D2","D4","D7")) 
saveTPM.bytime.tiss("Old", days =c("D0","D2","D4","D7"))

# calculate Tau for "meanTPMalldays
tau.bytimeandtype <- function( fileprefix,age,extension){
  t <- read.table(paste0(fileprefix,age, extension), 
                  sep='\t', header=T)
  print(head(t))
  rownames(t) <- t$id
  t$id <- NULL
  log.t <- log10(t+1)
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

