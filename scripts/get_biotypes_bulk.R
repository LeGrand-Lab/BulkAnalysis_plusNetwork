library(dplyr)
library(ggplot2)
library(tidyverse)
library(cowplot)
library(biomaRt)
library("AnnotationHub")
library(EnsDb.Mmusculus.v79)
###
# initial_expl.R has produced 'geneID_bc_matrix.rds'.
# here import it to get the biotype information
##
setwd("~/bulk_analysis/")
fileout <- "data/biotype_bulk.csv"
df <- readRDS("data/geneID_samp_matrix.rds")


# GET BIOTYPE INFORMATION
# =======================================================================
# WITH BIOMART
mmusRef <- useMart("ENSEMBL_MART_ENSEMBL", dataset="mmusculus_gene_ensembl",
                   host="www.ensembl.org")
infmart <- getBM(attributes = c("external_gene_name", "ensembl_gene_id",
                                "gene_biotype","strain_name", "chromosome_name"), 
                 values = rownames(df),
                 filters = "ensembl_gene_id",
                 mart=mmusRef)
#pop = infmart %>% group_by(strain_name) %>% summarise(count = n())
# CUSTOM DEDUPLICATE: split(uniq, dup), then exclude non canonical chr from dup:
uniqmart <- infmart %>% group_by(ensembl_gene_id) %>% dplyr::filter(n()==1) 
dupmart <- infmart %>% group_by(ensembl_gene_id) %>% 
  dplyr::filter(n()>1) %>%
  dplyr::filter(!str_detect(chromosome_name,'\\.')) %>%
  dplyr::filter(!str_detect(chromosome_name, "_")) 
# retain first occurence (sorted desc biotype)
dupmart <- dupmart %>% arrange(desc(gene_biotype)) %>%
  group_by(ensembl_gene_id) %>% 
  dplyr::filter(row_number(strain_name)==1)
furtherdup <- dupmart %>% group_by(ensembl_gene_id) %>% dplyr::filter(n()>1)
if (dim(furtherdup)[1]==0){ ("no more duplicates in dupmart") }

#new infmart (merge uniq and dup)
infmart <- rbind(uniqmart,dupmart)
infmart$source <- "MART_ENSEMBL"
infmart$strain_name[infmart$strain_name==""] <- "canonical" #fill if empty

## take strains datasets:
allstrains <- listDatasets(useMart(host = "www.ensembl.org", biomart = "ENSEMBL_MART_MOUSE"))$dataset
# give good format to allstrains vector (relocating dba as first, ...etc)
if (allstrains[1] != "mmdba2j_gene_ensembl"){
  allstrains[(match("mmdba2j_gene_ensembl",allstrains))] <- allstrains[1]
  allstrains[1] <- "mmdba2j_gene_ensembl"
  names(allstrains) <- c("DBA/2J",  "A/J", "AKR/J", "BALBC/J",    
                         "C3H/HEJ", "C57BL/6NJ", "CAST/EIJ", "CBA/J" ,     
                         "129S1/SVIMJ", "FVB/NJ", "LP/J", "NOD/SHILTJ", 
                         "NZO/HLLTJ", "PWK/PHJ","WSB/EIJ")
}
# for not yet found genes, create empty tibble to fillup with strains datasets 
absent = rownames(df)[!rownames(df) %in% infmart$ensembl_gene_id]
dfstrains <- infmart %>% dplyr::filter(n()>5) 
# parse all mart datasets:
for (i in 1:length(allstrains)){
  x = names(allstrains)[i]
  martx = useMart("ENSEMBL_MART_MOUSE", dataset=allstrains[x], host="www.ensembl.org")
  infx <- getBM(attributes = c( "external_gene_name","ensembl_gene_id",
                                "gene_biotype", "chromosome_name"), 
                values = absent, filters = "ensembl_gene_id",
                mart=martx)
  #fill strain (not by default in _MOUSE mart)
  if (dim(infx)[1]>0){
    infx <- infx %>% dplyr::filter(ensembl_gene_id %in% absent) %>%
      mutate(strain_name=x) %>% relocate(strain_name, .before="chromosome_name")
    infx <- infx %>% arrange(desc(gene_biotype)) %>% 
      group_by(ensembl_gene_id) %>% 
      dplyr::filter(row_number(strain_name)==1)  #deduplicate 
    for (i in colnames(infx)){ infx[[i]] <- as.character(infx[[i]]) }
    infx$source <- "ENSEMBL_MART_MOUSE"
    dfstrains <- dplyr::bind_rows(dfstrains,infx)
    tmp = absent[!absent %in% infx$ensembl_gene_id]
    absent = tmp  # update absent
  }
  print(paste("after",x, "dataset matching, absent left:" ))
  print(length(absent))
}


# WITH annotationdbi
# keytypes(EnsDb.Mmusculus.v79)
df_dbi <- select(EnsDb.Mmusculus.v79, keys=absent, keytype="GENEID",
                 columns=c("SYMBOL","GENEID","GENEBIOTYPE", "SEQNAME"))

df_dbi <- df_dbi %>% mutate(strain_name = "canonical") %>% 
  relocate(strain_name,.before=SEQNAME)
df_dbi$source = "AnnotationDBI"
colnames(df_dbi) <- colnames(dfstrains) # appropriate colnames

# NA's dataframe: those not found with any of automatized queries performed:
finalabsent <- absent[!absent %in% df_dbi$ensembl_gene_id]
t=length(finalabsent)
nas_df <- data.frame("external_gene_name"=finalabsent, "ensembl_gene_id"=rep("NOTFOUND",t),
                     "gene_biotype"=rep("NOTFOUND",t), "strain_name"=rep("NOTFOUND",t),
                     "chromosome_name"=rep("NOTFOUND",t), "source"=rep("NOTFOUND",t))


# merge the three dataframes and save into csv and SAVE:
finaldf <- dplyr::bind_rows(infmart, dfstrains,df_dbi,nas_df)
write.table(finaldf, file=fileout, sep=';', col.names = T, row.names = F)

# =====   END










