# prepare TPM matrix for natmi
library("dplyr")
library("tidyverse")

setwd("~/bulk_analysis/")
odir = "data/"
metadata.rds <- "data/metadata.rds"
prefil_tpm <- "data/prefiltered_TPM.rds"

fTPM <- readRDS(prefil_tpm)
metadata <- readRDS(metadata.rds)

TPM.y <- fTPM[, str_detect(colnames(fTPM),"Young")]
TPM.o <- fTPM[, str_detect(colnames(fTPM),"Old")]

# split into young and old, and filter out very low TPM:
keep.y <- apply(TPM.y, 1, function(row) sum(row >= 0.001) >= 3 )
keep.o <- apply(TPM.o, 1, function(row) sum(row >= 0.001) >= 3 )
TPM.y <- TPM.y[keep.y, ]
TPM.o <- TPM.o[keep.o, ]
dim(TPM.y)

## working only  'young'  for the moment:
# lets work on one single day first (D7), 
day = "D7"
age = "Young"
mymx <- TPM.y[,str_detect(colnames(TPM.y),day)]
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
  # calculate geometric mean across genes (row:
  draft_list[[k]] <- apply(tempo,1,function(row){prod(row)^(1/n.samples)})
}
newmx <- as.data.frame(draft_list)
#  filter out zeros TPM and save table
keepfinal <- apply(newmx[-1], 1, function(row) sum(row)>0)
newmx <- newmx[keepfinal,]
write.table(newmx, paste0(odir, paste("meanTPM",age,day,".txt"), sep='\t',
            col.names = T, row.names = F))

# loop to generalise to desired date , use young matrix first:
days <- c("D0","D2", "D4")
age="Young"
for (day in days){
  print(day)
  mymx <- TPM.y[,str_detect(colnames(TPM.y),day)]
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
    # calculate geometric mean across genes (row:
    draft_list[[k]] <- apply(tempo,1,function(row){prod(row)^(1/n.samples)})
  }
  newmx <- as.data.frame(draft_list)
  #  filter out zeros TPM and save table
  keepfinal <- apply(newmx[-1], 1, function(row) sum(row)>0)
  newmx <- newmx[keepfinal,]
  write.table(newmx, paste0(odir, paste0("meanTPM",age,day,".txt"), sep='\t',
                            col.names = T, row.names = F))
}

# now use old matrix
days <- c("D0","D2", "D4", "D7")
age="Old"
for (day in days){
  print(day)
  mymx <- TPM.o[,str_detect(colnames(TPM.o),day)]
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
    # calculate geometric mean across genes (row:
    draft_list[[k]] <- apply(tempo,1,function(row){prod(row)^(1/n.samples)})
  }
  newmx <- as.data.frame(draft_list)
  #  filter out zeros TPM and save table
  keepfinal <- apply(newmx[-1], 1, function(row) sum(row)>0)
  newmx <- newmx[keepfinal,]
  write.table(newmx, paste0(odir, paste0("meanTPM",age,day,".txt"), sep='\t',
                            col.names = T, row.names = F))
}
