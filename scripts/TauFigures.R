## for figure purposes
# TODO: improve to yield ggplots !!!!
calculateTau <- function(vec){
  if(max(vec) == 0){ # avoid zero division error
    tau.index <- 0
  }else{  
    tau.index <- (sum(1-(vec/max(vec))))/(length(vec)-1)
  }
  return(tau.index)
}

savee <- function(age, cts = list(), days=c("D0")){
  for (i in days){
    ktab <- read.table(paste0("data/meanTPM_",age,i,".txt"), sep='\t', header=T,
                       row.names=1) 
    logtab <- log10(ktab+1) 
    hist(unlist(logtab), col="lightblue")
    q25 = quantile(unlist(logtab), 0.25)
    keep <- apply(logtab, 1, function(x) sum(x >= 0) == length(x) &
                    sum(x > q25) > 1)  #and at least one over median
    logtab <- logtab[keep,] # for tau calc purposes
    hist(unlist(logtab), col="cyan", main="filterdone")
    print(dim(logtab))
    tau_res <- tibble("id"=rownames(logtab))
    tau_res$symbol <- genes_df[match(rownames(logtab),genes_df$Geneid),]$symbol
    tau_res$Tau <- apply(logtab, 1, function(row) calculateTau(row))
    hist(tau_res$Tau, col="wheat", main=paste(age,": ",i))
  }
  return("tau calculated & txt files saved into 'Tau/'")
}
print( savee("Young", days=c("D0","D2","D4","D7")) )
print( savee("Old",  days=c("D0","D2","D4","D7")) ) 

summary(unlist(logtab))
