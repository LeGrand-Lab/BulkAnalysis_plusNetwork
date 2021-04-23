# Prepares TPM matrices for :
# - Natmi Ligand Receptor network 

library("dplyr")
library("tidyverse")

setwd("~/BulkAnalysis_plusNetwork/")
odir = "inDataNatmi/"
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

ages <- c("Young","Old")
days <- c("D0","D2", "D4", "D7")
for (iage in ages){
  TPM.i <- fTPM[, str_detect(colnames(fTPM),iage)]
  meta.i <- metadata %>% filter(age==iage)
  for (day in days){
    print(day)
    premx <- TPM.i[,str_detect(colnames(TPM.i), day)]
    mymx <- tibble("gene"=rownames(premx))
    mymx <- cbind(mymx,premx)
    metax <- meta.i %>% filter(time==day)
    print(all(colnames(mymx)[-1]==rownames(metax)))
    annot <- metax %>% mutate(Sample=newname, Cell_type=type) %>%
      select(Sample,Cell_type) %>% tibble()
    write.table(mymx, paste0(odir,"TPM_",iage,day,".txt"), sep='\t', 
              col.names = T, row.names = F)
    write.table(annot, paste0(odir,"annot_",iage,day,".txt"),sep='\t', 
              col.names = T, row.names = F)
  }
}

#### END