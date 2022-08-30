###
# Takes multiple files from StringTie results
# outputs expression csv and counts rds files
# -
# johaGL 2021
##

#library(reshape2)
library(tidyverse)
library(stringr)

setwd("~/BulkAnalysis_plusNetwork/")
data_p <- "/media/bioinfo/Stockage/DATA_RNAageing/DATA/"
m_p <- "merged_gene_counts/" # /media/bioinfo/Stockage/DATA_RNAageing/total_results_nfcore/featureCounts/gene_counts
d_p <- "StringtieFPKM_gene_abund_all_batch/"

odir = "data/"
#outputs:
g_samp_mat = "geneID_samp_matrix.rds"
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
    tmpdf <- dplyr::select(tmpdf,-Number)
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

write.table(batchesinfo, paste0(odir,batch_f), col.names = T, sep="\t")
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
write.table(genesinfos, paste0(odir,genes_f), col.names=T,sep="\t")

# =========================================================================
# transform data.frame into matrix
# =========================================================================
rownames(df) <- df$Geneid
df <- dplyr::select(df,-c(Geneid,gene_name))
saveRDS(df, paste0(odir,g_samp_mat))
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



# extract only genes infos first into 'featdf': 
fullpath = paste0(data_p, d_p, tie_fl[1])
tmp <- read.table(fullpath, sep="\t", header=T)
featdf <- tmp %>% dplyr::select(!c(Coverage,FPKM,TPM))
for (i in 1:length(tie_fl)){
  print(tie_fl[i])
  fullpath = paste0(data_p, d_p, tie_fl[i])
  pre_feat <- tmp %>% dplyr::select(!c(Coverage,FPKM,TPM))
  featdf <- rbind(featdf, pre_feat)
  featdf <- unique(featdf) #remove duplicated rows
}
#verify uniqueness in this stringtie Gene.ID column:
print(all(featdf$Gene.ID==unique(featdf$Gene.ID)))  # FALSE 
# Gene.ID is NOT unique
noccurs = data.frame(table(featdf$Gene.ID))
noccurs[noccurs$Freq >1,]   
# Var1 Freq
# 5092  ENSMUSG00000025515    5
# 5178  ENSMUSG00000025762    2
# 5388  ENSMUSG00000026162    2
# 16072 ENSMUSG00000056050    2
# 32915 ENSMUSG00000093803    4
# 38813 ENSMUSG00000100826    4

# Each oune of these genes has start and end positions different by some bases:
# Gene.ID Gene.Name Reference Strand     Start       End
# 2184 ENSMUSG00000056050      Mia3         1      - 183326725 183354040
# 2185 ENSMUSG00000056050      Mia3         1      - 183356346 183369553
write.table(featdf, paste0(odir,outfeat), sep="\t", col.names = T, row.names = F)
rowsfull = unique(featdf$Gene.ID)
# for simplicity, took unique Gene.ID from this df

TPM = data.frame("Gene.ID"=rowsfull)
FPKM = data.frame("Gene.ID"=rowsfull)
cov = data.frame("Gene.ID"=rowsfull)

rownames(cov) = rowsfull
rownames(TPM) = rowsfull
rownames(FPKM) = rowsfull

for (i in 1:length(tie_fl)){
  print(tie_fl[i])
  fullpath = paste0(data_p, d_p, tie_fl[i])
  tmp <- read.table(fullpath, sep="\t", header=T)
  # TPM
  preTPM <- tmp$TPM
  names(preTPM) <- tmp$Gene.ID
  TPM[[newnames[i]]] = preTPM[match(rownames(TPM),names(preTPM))]
  # fpkm
  preFPKM <- tmp$FPKM
  names(preFPKM) <- tmp$Gene.ID
  FPKM[[newnames[i]]] = preFPKM[match(rownames(FPKM),names(preFPKM))]
  #coverage:
  pre_cov <- tmp$Coverage
  names(pre_cov) <- tmp$Gene.ID
  cov[[newnames[i]]] = pre_cov[match(rownames(cov),names(pre_cov))]
}
TPM[is.na(TPM)] <- 0
FPKM[is.na(FPKM)] <- 0
cov[is.na(cov)] <- 0

#check ok colnames:
print(all(colnames(TPM)[-1] == newnames))

write.table(TPM, paste0(odir,outTPM), sep="\t", col.names = T, row.names = F)
write.table(FPKM, paste0(odir,outFPKM), sep="\t", col.names = T, row.names = F)
write.table(cov, paste0(odir,outcov), sep="\t", col.names = T, row.names = F)

#outTPM 
#outFPKM 
#outcov
#outfeat

# cool hint: subset only rows containing NA
# df[!complete.cases(df),] 