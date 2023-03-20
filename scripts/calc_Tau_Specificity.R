##
# Tau specificity index calculation
# --
# johaGL 2021
#
# This script will save into different txt and xlsx files to 'Tau/' folder.
# in this way, separately for young and old:
#       ids    symbol    Tau   class  whichMAX   nbMAX    day
# if several tissues exhibit max logTPM, they are all registered separated 
# by commas as single string in column whichMAX
# so 'housekeeping' and 'specific' are going to appear at the  class column
# Tau >= 0.8 ~ "specific",
# Tau >= 0.5 & Tau < 0.8 ~ "intermediate",
# Tau < 0.5 ~ "housekeeping"
# A consensus is also generated : 
library(dplyr)
library(tidyverse)
library(openxlsx)
library(reshape2)

setwd("~/BulkAnalysis_plusNetwork2/")

genes_df <- read.table("data/genesinfo.csv",sep="\t",header=T)
consensusfile <- "conseTau_ensemblid.rds"
needconsensus <- T
###  CALCULATE
# we want, for type&day,´ to obtain sample specific genes, 
# and done separately by age (old separated from young)
# import TPM matrices by day by age
# do Tau by cell type but also by celltype&day
# ====================== Calculate =============================
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
saveTau.bytissue <- function(age,  days, myquantile=0.25){
  for (i in days){
    ktab <- read.table(paste0("data/meanTPM_",age,i,".txt"), sep='\t', header=T,
                       row.names=1) 
    logtab <- log10(ktab+1) 
    # get rid homogeneously low rows in matrix, setting a min quantile cutoff
    q.cutoff <- quantile(unlist(logtab),myquantile)
    keep <- apply(logtab, 1, function(x) sum(x >= 0) == length(x) &
                    sum(x > q.cutoff) >= 1)  #and at least one over this cutoff
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
  return("tau calculated ==> txt files saved into 'Tau/'")
}
print( saveTau.bytissue("Young", days=c("D0","D2","D4","D7")) )
print( saveTau.bytissue("Old",  days=c("D0","D2","D4","D7")) ) 

## ===================== Filter generated Tau matrices =======================
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
    #write.table(itab, paste0("Tau/specifivOnly/TauSpecificity_",age,day,".txt"), sep='\t',
    #             col.names=T, row.names=T)
    # prepare xls file
    xfile = paste0("Tau/specificOnly/TauSpecificity_",age,day,".xlsx")
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


# END formal Tau calc
## ============================ Consensus ====================================

if (needconsensus){ 
  print("building consensus")
  consensus_tau <- list()
  daysv = c("D0","D2","D4","D7")
  for (day in daysv){
    ty <- read.table(paste0("Tau/TauSpecificity_Young",day, ".txt"), sep="\t", header=T)
    to <- read.table(paste0("Tau/TauSpecificity_Old",day, ".txt"), sep="\t", header=T)
    # remember = Tau is calculated on meanTPM (mean over replicates' TPMs)
    #  therefore, its value is also influenced by changes in expression due to age
    # for example, a hypothetical gene1 : [FAPS, sCs, ECs, M1], their meanTPMs
    # being :   Old =  [100,100,100,100]  ; Young = [100,2000,100,100]
    # check documentation for justification of following code: 
    ty$age = "Young"
    to$age = "Old"
    rownames(ty) = ty$id
    rownames(to) = to$id
    tmpdf = data.frame( "id"=character(), "symbol"=character(), "Tau"=numeric(), 
                        "class"=character(),  "whichMAX"=character(), 
                        "nbMAX"=numeric(), "exclusiveOld"=numeric() )
    allids.here = unique(c(ty$id,to$id))
    for (i in allids.here){
      if (is.na(ty[i,]$id)){
        print(paste0("detected exclusive old : ", i))
        tmpdf[i,] <- c( i, to[i,]$symbol,  to[i,]$Tau, 
                        to[i,]$class, to[i,]$whichMAX, to[i,]$nbMAX, 1  )
      } else {
        tmpdf[i, ] <-  c( i, ty[i,]$symbol, ty[i,]$Tau, 
                          ty[i,]$class, ty[i,]$whichMAX, ty[i,]$nbMAX, 0 )
      }   
    }
    consensus_tau[[day]] = tmpdf
  }
  saveRDS(consensus_tau, paste0("Tau/", consensusfile))
}else{
  print(paste0("consensuslist already exists in :",
               "Tau/", consensusfile)) }

# ================================ * * ====================================
# =============================== ADDENDUM ================================
# =========================================================================
# I. PROCESSING both by day and by age TPM matrices for eventual
# downstream purposes :
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
  write.table(fussi, paste0("data/future4Tau/meanTPMalldays_",age,".txt"), 
              sep='\t', col.names=T, row.names=T)
}

saveTPM.bytime.tiss("Young", days=c("D0","D2","D4","D7")) 
saveTPM.bytime.tiss("Old", days =c("D0","D2","D4","D7"))

# END processing eventual interesting datasets

# II.  EXPERIMENTAL : # calculate Tau for "meanTPMalldays.." generated by I.
tau.bytimeandtype <- function( fileprefix,age,extension, myquantile=0.25){
  t <- read.table(paste0(fileprefix,age, extension), sep='\t', header=T)
  rownames(t) <- t$id
  t$id <- NULL
  log.t <- log10(t+1)
  # get rid homogeneously low rows in matrix, setting a min quantile cutoff
  q.cutoff <- quantile(unlist(log.t),myquantile)
  keep <- apply(log.t, 1, function(x) sum(x >= 0) == length(x) &
                  sum(x > q.cutoff) >= 1)  #and at least one over this value
  log.t <- log.t[keep,] 
  tissues <- colnames(log.t)
  tau_res <- tibble("id"=rownames(log.t))
  tau_res$symbol <- genes_df[match(rownames(log.t),genes_df$Geneid),]$symbol
  tau_res$Tau <- apply(log.t, 1, function(row) calculateTau(row))
  tau_res <- tau_res %>% mutate(class = case_when(
    Tau >= 0.8 ~ "specific",
    Tau >= 0.5 & Tau < 0.8 ~ "intermediate",
    Tau < 0.5 ~ "housekeeping"
  ))
  whichMAX <- c()
  nbMAX <- c()
  maxlog10TPM <- c()
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
    maxlog10TPM <- c(maxlog10TPM, max(log.t[j,]) ) 
  }
  tau_res <- cbind(tau_res, whichMAX,nbMAX, maxlog10TPM)
  write.table(tau_res, paste0("Tau/Tau_TimeType/TauSpecificity_TimeType_",
                     age,".txt"), 
              sep='\t', col.names=T, row.names=T)
  return("ok")
}
age = "Old"
tau.bytimeandtype( "data/future4Tau/meanTPMalldays_", age,".txt")
age= "Young"
tau.bytimeandtype( "data/future4Tau/meanTPMalldays_", age,".txt")

# =========================================================================
## END ADDENDUM