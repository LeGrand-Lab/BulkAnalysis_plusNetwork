# prepare TPM matrices for natmi
# by age , by day

library("dplyr")
library("tidyverse")

setwd("~/BulkAnalysis_plusNetwork/")
odir = "data/"
metadata.rds <- "data/metadata.rds"
prefil_tpm <- "data/prefiltered_TPM.rds"

fTPM <- readRDS(prefil_tpm)
metadata <- readRDS(metadata.rds)

# check if coverage is ok :
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

# filter out zero coverage rows from fTPM
#  see how many are null coverage in fTPM
badgenes <- rownames(fTPM)[rownames(fTPM) %in% rownames(genesZeroCoverage)]
length(badgenes)
fTPM <- fTPM[!rownames(fTPM) %in% idszerocover,] 

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
      samples <- names(cellty[cellty==k])
      n.samples <- length(samples) 
      tempo <- as_tibble(mymx) %>% select(samples) 
      # calculate geometric mean across genes (row)
      draft_list[[k]] <- apply(tempo,1,function(row){prod(row)^(1/n.samples)})
    }
    newmx <- as.data.frame(draft_list)
    #  filter out zeros TPM and save table
    keepfinal <- apply(newmx[-1], 1, function(row) sum(row)>0)
    newmx <- newmx[keepfinal,]
    write.table(newmx, paste0(odir,"meanTPM_",age,day,".txt"), sep='\t', 
                col.names = T, row.names = F)
  }
}




## first draft:
# ## working only  'young'  for the moment:
# TPM.y <- fTPM[, str_detect(colnames(fTPM),"Young")]
# # lets work on one single day first (D7), 
# day = "D7"
# age = "Young"
# mymx <- TPM.y[,str_detect(colnames(TPM.y),day)]
# # calculate mean TPM by celltype : get celltypes present in this data
# # by using metadata
# cellty <- metadata[match(colnames(mymx),metadata$newname),]$type
# names(cellty) <- colnames(mymx)
# 
# ## loop across celltypes
# uniqcellty <- unique(cellty)
# draft_list <- list()
# draft_list[["ensembl"]] <- rownames(mymx)
# for (k in uniqcellty){
#   samples <- names(cellty[cellty==k])
#   n.samples <- length(samples) 
#   tempo <- as_tibble(mymx) %>% select(samples) 
#   # calculate geometric mean across genes (row:
#   draft_list[[k]] <- apply(tempo,1,function(row){prod(row)^(1/n.samples)})
# }
# newmx <- as.data.frame(draft_list)
# #  filter out zeros TPM and save table
# keepfinal <- apply(newmx[-1], 1, function(row) sum(row)>0)
# newmx <- newmx[keepfinal,]
# write.table(newmx, paste0(odir, paste("meanTPM",age,day,".txt"), sep='\t',
#                           col.names = T, row.names = F))