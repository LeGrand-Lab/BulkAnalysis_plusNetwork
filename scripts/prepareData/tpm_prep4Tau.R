# Prepares TPM matrices for * Tau * specificity index
# As biological replicates will be subjected to geometric mean, CAUTION ! :
# when using **geometric mean** : 
# MAKE SURE NO VALUES == ZERO ARE PRESENT IN THE MATRIX
# JohaGL 2021
# Inspired From M. Love (estimateSizeFactors{DESeq2}), about geomMean : 
 # """sample with a zero (a problem for the default method,
 # as the geometric mean becomes zero, and the ratio undefined)"""

library("dplyr")
library("tidyverse")

setwd("~/BulkAnalysis_plusNetwork/")
odir = "data/"
metadata.rds <- "data/metadata.rds"
prefil_tpm <- "data/prefiltered_TPM.rds"

fTPM <- readRDS(prefil_tpm)
metadata <- readRDS(metadata.rds)

# check if coverage is ok :
if (!file.exists(paste0(odir,"nullCoverageEnsemblids.txt"))){
  print("find zero coverage genes, save into data/")
  coverage <- read.table("data/COVmatrix.csv", sep="\t", header=T)
  rownames(coverage) <- coverage$Gene.ID
  coverage$Gene.ID <- NULL
  covCOLNAMES <- metadata[match(colnames(coverage), metadata$sample),]$newname
  coverage <- coverage[,!is.na(covCOLNAMES)] # drop out bad samples
  max(as.numeric(unlist(coverage)))  #[1] 270906.8
  genesZeroCoverage <- coverage %>% filter(rowSums(coverage) == 0)
  # save null coverage (i mean, ZERO across ALL samples) ensembl ids vector :
  idszerocover <- rownames(genesZeroCoverage)
  write_lines(idszerocover, 
              file=paste0(odir,"nullCoverageEnsemblids.txt"), sep='\n')
}
idszerocover = read_lines(file=paste0(odir,"nullCoverageEnsemblids.txt"))
# filter out zero coverage rows from fTPM
fTPM <- fTPM[!rownames(fTPM) %in% idszerocover,] 

## explore means by gene
# =========================================================================
# take one gene at random
#  samptpm <- fTPM[sample(nrow(fTPM), size=1),] # row at random
# mediu <- apply(fTPM, 1, function(x) sum(x)> 2000 & sum(x) < 5000)
# > sample(rownames(fTPM[mediu,]), 1)
kgenes <- c("ENSMUSG00000058794", "ENSMUSG00000028059","ENSMUSG00000002944")
kplots <- list()
#highly expressed, medium and low expressed genes, respectively
for (kgene in kgenes){
  One <- metadata %>% dplyr::select(newname, time, age, type)
  One$value <- fTPM[kgene,]
  One <- One %>% mutate(mycolor = ifelse(age=="Old","black","gray"))
  repli <- sapply(One$newname, function(x) str_split(x, "_")[[1]][2])
  categ <- sapply(One$newname, function(x) str_split(x, "_")[[1]][1])
  One$categ <- categ
  prep <- data.frame("categ"=unique(categ),
                     geommean = NA,
                     arithmean =NA,
                     median = NA)
  rownames(prep) <- prep$categ
  for (i in rownames(prep)){
    thisvals <- One %>% filter(categ==prep[i,"categ"]) %>% pull(value)
    nozerovals <- thisvals
    nozerovals[nozerovals==0] <- 0.00001  # initial test with this ZERO treatment 
    prep[i,"geommean"] = prod(nozerovals)^(1/length(nozerovals))
    prep[i, "arithmean"] = mean(thisvals)
    prep[i, "median"] = unname(quantile(thisvals, 0.5))
  }
  prep2 <- inner_join(prep, One, by="categ")
  g <- ggplot(One, aes(x=age, y=value, fill=age)) +
    geom_dotplot(binaxis='y', dotsize = 3, color="white", alpha=.3,
                 show.legend=FALSE   ) +
    geom_point(shape=18, data=prep2,aes( x=age, y=arithmean), 
               color="purple" , size= 3, alpha=.3)+
    geom_point(shape=18, data=prep2,aes( x=age, y=geommean), 
               color="black", size = 3, alpha=.3 )+
    geom_point(shape=9, data=prep2,aes( x=age, y=median), 
               color="darkblue",size = 2, alpha=.3)+
    facet_grid(vars(time), vars(type)) +
    labs(title="Testing geometric mean", subtitle=kgene,
         caption= "crossedDiamond = median; purpleDiamond = arithmethic mean; 
       blackDiamond = geometric mean.   *Note: ignore righthand legend") +
    theme(legend.position="none") +
    theme_bw()
  kplots[[kgene]] <- g
}
pdf(paste0("plotsPrelim/tpmmetrics_randomGenes.pdf"), width=8)
kplots[[1]]
kplots[[2]]
kplots[[3]]
dev.off()
# ============================================================================
# end explore means by gene

## DISCOVER MIN VAL FOR REPLACING ZEROES
print("exploring minimal existing value not being ZERO TPM:")
hasNOzeros = apply(fTPM, 1, function(row) sum(row > 0) ==length(row))
lowestTPM = min(unlist(fTPM[hasNOzeros,]))
print(paste("this is the min value admited, will replace zeroes",lowestTPM))
print('replacing zero values by lowest found in entire matrix')
fTPM[fTPM < lowestTPM] <- lowestTPM


ages <- c("Young","Old")
days <- c("D0","D2", "D4", "D7")
for (age in ages){
  TPM.i <- fTPM[, str_detect(colnames(fTPM),age)]
  for (day in days){
    print(day)
    mymx <- TPM.i[,str_detect(colnames(TPM.i), day)]
    keep.x <- apply(mymx, 1, function(row) sum(row >= 0.001) >= 3 )
    mymx <- mymx[keep.x, ]
    # calculate mean TPM by celltype : get celltypes present in this data
    # by using metadata
    cellty <- metadata[match(colnames(mymx),metadata$newname),]$type
    names(cellty) <- colnames(mymx)
    ## loop across celltypes
    uniqcellty <- unique(cellty)
    draft_list <- list()
    draft_list[["ensembl"]] <- rownames(mymx)
    for (k in uniqcellty){
      biolreplic <- names(cellty[cellty==k]) 
      n.bioreplic <- length(biolreplic) 
      colsbiorep <- mymx[,biolreplic]
      # calculate geometric mean across genes (row = biol replicates TPM values)
      draft_list[[k]] <- apply(colsbiorep,1,function(row){
        prod(row)^(1/n.samples)
        })
    }
    newmx <- as.data.frame(draft_list)
    #  filter out zeros TPM and save table
    keepfinal <- apply(newmx[-1], 1, function(row) sum(row)>0)
    newmx <- newmx[keepfinal,]
    write.table(newmx, paste0(odir,"meanTPM_",age,day,".txt"), sep='\t', 
                col.names = T, row.names = F)
  }
}

#### END