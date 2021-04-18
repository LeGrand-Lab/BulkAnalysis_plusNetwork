##  closely related to 'calc_Tau_Specificity.R'
# this script is only for figure purposes

library(dplyr)
library(tidyverse)
library(openxlsx)
library(reshape2)

setwd("~/BulkAnalysis_plusNetwork/")


calculateTau <- function(vec){
  if(max(vec) == 0){ # avoid zero division error
    tau.index <- 0
  }else{  
    tau.index <- (sum(1-(vec/max(vec))))/(length(vec)-1)
  }
  return(tau.index)
}


age =  "Young"
#age = "Old"
#see tau plots
pdf(paste0("Tau/","ZplotDetailscalc_", age,".pdf"))
days=c("D0","D2","D4","D7")
par(mfrow = c(4,3) )
  for (i in days){
    ktab <- read.table(paste0("data/meanTPM_",age,i,".txt"), sep='\t', header=T,
                       row.names=1) 
    logtab <- log10(ktab+1) 
    hist(unlist(logtab),xlim=c(-3,8), col="wheat", prob=T, cex.main=.9,
         main=paste("compare logarithmic transform\n (TPM+1)",i, age), 
         xlab = "bars:log10 (lines:log2) TPM+1")
    lines(density(unlist(log(ktab+1))),  col="red", cex.main=.7)
    text(4,0.15,"log base 2", col="red")
    # get rid homogeneously low rows in matrix, setting a min quantile cutoff
    myquantile = 0.25
    q.cutoff <- quantile(unlist(logtab),myquantile)
    keep <- apply(logtab, 1, function(x) sum(x >= 0) == length(x) &
                    sum(x > q.cutoff) >= 1)  #and at least one over this value
    logtab <- logtab[keep,] 
    hist(unlist(logtab), col="dodgerblue", 
         main=paste("Filtered rows log10(TPM+1)\n",
        "having no values over quantile", myquantile, "(",round(q.cutoff,4),")"  ),
         cex.main=.7, xlab = "log10(TPM+1)")
    print(dim(logtab))
    tau_res <- tibble("id"=rownames(logtab))
    tau_res$Tau <- apply(logtab, 1, function(row) calculateTau(row))
    hist(tau_res$Tau, col="orange", main=paste("Tau index at", i,",", age),
         xlab= 'Tau index')
    
  }
par(mfrow=c(1,1))
dev.off()


