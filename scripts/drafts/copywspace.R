library(dplyr)
library(tidyverse)
library(Matrix)
setwd("~/bulk_analysis/")

#design has been copied into this wdir:
data_p <- "/media/bioinfo/Stockage/DATA_RNAageing/DATA/"
design <- read.table("design.csv",header=T,sep=",", row.names=1)

mat <- readRDS("geneID_bc_matrix.rds")

absentbc <- setdiff(colnames(mat),design$sample)
# TODO: why are these samples not included in design?????:
#[1] "X2.422363" "X2.422364" "X2.424030" "X2.424031" "X2.428151" "X2.428100" "X2.429794" "X2.429793"
# maybe just "negative" controls, lets filter matrix to check if thse samples remain


# =========================================================================
# TODO
# =========================================================================
# 1. import TPM information:
#    -calculate mean TPM, filter out genes < 1 mean TPM
TPM <- read.table("TPMmatrix.csv", header=T,sep="\t",row.names=1)
head(TPM)
meanTPM <- apply(TPM,1,mean) # or sumsTPM/dim(TPM)[2]
ggtpm1 <- ggplot(as_tibble(meanTPM),aes(log10(meanTPM+1e-15))) + 
  ggtitle("mean TPM distribution across genes") + 
  geom_histogram(aes(y=..density..), color="white",fill="#dd9c00",alpha=.7) + 
  geom_vline(xintercept=log(1),linetype = "dashed" ) +
  geom_density( color="gold", alpha=.2) +
  theme_bw()

# we calculated mean TPM per gene across all samples, here we call it 'TPM value':
# the distribution of the log for these values indicates that aprroximately:
# 20 % of values correspond to  Zero TPM 
# ~80% of TPM values are comprised beetween 12 and 0.02 
# most frequent TPM value is 1
# so NO! cutoff =1 for us here.
# LETS SIMULATE CANDIDATE THRESHOLDS (logarithmic) and do plot
# Average TPM rule:
logtpmcutoff <- seq(from=-16, to=6, length=100)
cutoff_ngenes <- tibble("logTPM_cutoff" = logtpmcutoff,
                        "how_many_genes" = sapply(logtpmcutoff, function(x){
                          subv <- meanTPM[meanTPM >= 10^x]
                          return(length(subv))  }) #end sapply
                      )

ggtpm2 <- ggplot(cutoff_ngenes) +
   geom_line(aes(x = logTPM_cutoff, y=how_many_genes)) +
  ggtitle("logTPM threshold vs genes preserved") +
  theme_bw()
# at least one half rule:
# how many genes have >= half TPMs values over threshold
cutoff_ngenes$at_least_half <- sapply(logtpmcutoff, function(x){
  sum(apply(TPM > x^10, 1 , mean) >= 0.5 )
})

test <- c()
for (i in 1:100){
  x = logtpmcutoff[i]^10
  tmp = c()
  for (j in 1:dim(TPM)[1]){
    avec <- TPM[j,][TPM[j,] > x]
    if (length(avec) >= 127/2){
      tmp <- c(tmp, TRUE)
    }else{tmp <- c(tmp, FALSE)}
  }
  test <- c(test,sum(tmp))
}

cutof_ngenes$bis_half <- test

ggplot(cutoff_ngenes) +
  geom_line(aes(x=logTPM_cutoff,y=bis_half)) +
  ggtitle("manual test") +
  theme_bw()


log_sel <- logtpmcutoff %>% filter(how_many_genes > 3000) %>% 
		desc(logTPM_cutoff) %>% top_n(1)

newtpm <- meanTPM[meanTPM>cutofftpm]
length(meanTPM[meanTPM>=0.08 & meanTPM<12])
ggtpm3 <- ggplot(as_tibble(),aes(log10(meanTPM))) + 
  ggtitle(paste0("mean TPM after genes filter, cutoff=",cutofftpm)) + 
  geom_histogram(aes(y=..density..), color="white",fill="#dd9c00",alpha=.7) + 
  geom_density( color="gold", alpha=.2) +
  theme_bw()

plot_grid(ggtpm1,ggtmp2, ggtpm3, ggtpm4, nrows = 2)


# 2. import information about biotype (done on matrix, see 'get_biotypes_bulk.R')
btf <- "biotype_bulk.csv"
biotype <- read.table(btf, sep=";",header=T)

# remove pseudogenes and ribo-related elements



# 3. remove rows of the DESeqDataSet that have no counts, 
# or only a single count across all samples.

# ###############################################################"
# # very first time exploration of samples
# names <- c()
# for (i in test){
#   tmpx <- read.table(paste0(data_p,m_p,i), sep='\t', header=T)
#   if(str_detect(i, "B2")){
#     # B2 : column 'Number' before Geneid and gene_name
#     prefixes <- sapply(colnames(tmpx)[c(-1,-2, -3)], function (y){
#       tmpsl <- str_split(y,"_")
#       return(tmpsl[[1]][1])}
#     )
#     names <- c(names,unname(prefixes))
#   }else{
#     prefixes <- sapply(colnames(tmpx)[c(-1,-2)], function (y){
#       tmpsl <- str_split(y,"_")
#       return(tmpsl[[1]][1])}
#     )
#     names <- c(names,unname(prefixes))
#   }
# }

