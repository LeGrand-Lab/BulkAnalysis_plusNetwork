library(reshape2)
library(tidyverse)
library(forcats) # tools categorical vars
library(UpSetR) #visuals intersecting sets
library(ggpubr)
library(kml)
library(ggplotify)
library(stringr)
library("DESeq2")
# library(umap)
setwd("~/bulk_analysis/")
data_p <- "/media/bioinfo/Stockage/DATA_RNAageing/DATA/"
m_p <- "merged_gene_counts/"
d_p <- "StringtieFPKM_gene_abund_all_batch/"


#outputs:
g_bc_mat = "geneID_bc_matrix.rds"
batch_f = "batchesinfo.csv"
genes_f = "genesinfo.csv"
outTPM <- "TPMmatrix.csv"
outFPKM <- "FPKMmatrix.csv"
outcov <- "COVmatrix.csv"
outfeat <- "genes_as_stringtie.csv"

  
# =========================================================================
# build dataframe from independent txt counts files
# =========================================================================
files <- list.files(paste0(data_p,m_p), pattern=".txt")
batchesinfo <- data.frame("sample"=character(),"batch"=character())

list_df <- list(length=length(files))
for(i in 1:length(files)){
  x <- files[i]
  tmpdf <- read.table(paste0(data_p, m_p, x), sep="\t", header=T)
  if (!str_detect(x,"B2")){
    origcols <- colnames(tmpdf)
    newcols <- sapply(origcols,function(coln)
      if(str_detect(coln,"merge")){
        str_split(coln,"_")[[1]][1]
        }else {return(coln)})
    colnames(tmpdf) <- unname(newcols)
    #print(head(tmpdf,3))
  }else{
    # B2: exceptional column 'Number' before Geneid and gene_name
    tmpdf <- select(tmpdf,-Number)
    # no need to format colnames, are already ok
    #print(head(tmpdf,3))
  }
  binfo <- data.frame(
    "sample"= colnames(tmpdf)[!colnames(tmpdf) %in% c("Geneid","gene_name")],
    "batch" = str_replace(str_split(x,"_")[[1]][4],".txt","")
  )
  batchesinfo <- rbind(batchesinfo,binfo)
  list_df[[i]] <- tmpdf
}

write.table(batchesinfo, batch_f, col.names = T, sep="\t")
df <- list_df[[1]]
df <- full_join(df,list_df[[2]], by=c("Geneid","gene_name"))
df <- full_join(df,list_df[[3]], by=c("Geneid","gene_name"))
df <- full_join(df,list_df[[4]], by=c("Geneid","gene_name"))

# verify sample name uniqueness
print(all(unique(batchesinfo$sample)==batchesinfo$sample))
# verify if any info leaks:
indepcols = 2 # the two first cols universal
for (k in 1:length(list_df)){
  indepcols = indepcols + (dim(list_df[[k]])[2])-2 #because id and name again
}
(indepcols)
(dim(df)[2])
rm(list_df)

# verify Geneid uniqueness and save gene_name information 
print(all(unique(df$Geneid)==df$Geneid))
genesinfos <- data.frame("Geneid"=df$Geneid, "symbol"=df$gene_name)
write.table(genesinfos, genes_f, col.names=T,sep="\t")

# =========================================================================
# transform data.frame into matrix
# =========================================================================
rownames(df) <- df$Geneid
df <- select(df,-c(Geneid,gene_name))
saveRDS(df, g_bc_mat)
rm(mat)

# =========================================================================
# GET: * TPM FPKM and Cov matrices *  and feature info dataframe
# =========================================================================

head(list.files(paste0(data_p,d_p)))
length(list.files(paste0(data_p,d_p)))
tie_fl <- list.files(paste0(data_p,d_p))
fullpath = paste0(data_p, d_p, tie_fl[1])
tmp <- read.table(fullpath, sep="\t", header=T)
# current sample format: "2-429798_mergeAligned.sortedByCoord.out.gene_abund.txt
# desired format: "X2.421544"
newnames <- unname(sapply(tie_fl, function(x) {
              newx <- str_split(x,"_")[[1]][1]
              newx <- paste0("X",str_replace(newx, "-", "."))}))
print(head(newnames))

cov <- tibble("Gene.ID"=tmp$Gene.ID)
cov[[newnames[1]]] = tmp$Coverage
TPM <- tibble("Gene.ID"=tmp$Gene.ID)
TPM[[newnames[1]]] = tmp$TPM
FPKM <- tibble("Gene.ID"=tmp$Gene.ID)
FPKM[[newnames[[1]]]] = tmp$FPKM
featdf <- tmp %>% select(!c(Coverage,FPKM,TPM))

for (i in 2:length(tie_fl)){
  print(tie_fl[i])
  fullpath = paste0(data_p, d_p, tie_fl[i])
  tmp <- read.table(fullpath, sep="\t", header=T)
  # TPM
  preTPM <- tibble("Gene.ID"=tmp$Gene.ID)
  preTPM[[newnames[i]]] = tmp$TPM
  TPM <- merge(TPM, preTPM, by="Gene.ID", all=T)
  # fpkm
  preFPKM <- tibble("Gene.ID"=tmp$Gene.ID)
  preFPKM[[newnames[[i]]]] 
  FPKM <- merge(FPKM, preFPKM, by="Gene.ID", all=T)
  #coverage:
  pre_cov <- tibble("Gene.ID"=tmp$Gene.ID)
  pre_cov[[newnames[i]]] = tmp$Coverage
  cov <- merge(cov, pre_cov, by="Gene.ID", all=TRUE) #all TRUE !! yes
  #feat
  pre_feat <- tmp %>% select(!c(Coverage,FPKM,TPM))
  featdf <- rbind(featdf, pre_feat)
  featdf <- unique(featdf) #remove duplicated rows
  print(dim(TPM))
}
TPM[is.na(TPM)] <- 0
FPKM[is.na(FPKM)] <- 0
cov[is.na(cov)] <- 0

#check ok colnames:
print(all(colnames(TPM)[-1] == newnames))

write.table(TPM, outTPM, sep="\t", col.names = T, row.names = F)
write.table(FPKM, outFPKM, sep="\t", col.names = T, row.names = F)
write.table(cov, outcov, sep="\t", col.names = T, row.names = F)
write.table(featdf, outfeat, sep="\t", col.names = T, row.names = F)
#outTPM 
#outFPKM 
#outcov
#outfeat

# cool hint: subset only rows containing NA
# df[!complete.cases(df),] 