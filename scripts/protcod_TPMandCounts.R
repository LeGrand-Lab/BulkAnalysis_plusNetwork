library(dplyr)
library(tidyverse)
library(cowplot)
library(forcats)
library(RColorBrewer)
library(MASS)
setwd("~/bulk_analysis/")

btf <- "data/biotype_bulk.csv"
mat <- readRDS("geneID_samp_matrix.rds")

#outputs
out.counts <- "data/protcod_counts.rds"
out.tpm <- "data/protcod_TPM.rds"
out1.pdf <- "plotsPrelim/pre_filtering_figures.pdf"
#rdatafile <- "onlyprot_TPMinf.RData"

# =========================================================================
# import matrices and visualize properties pre-filtering
# =========================================================================
# 1. import TPM matrix and raw count matrix:
TPM <- read.table("data/TPMmatrix.csv", header=T,sep="\t",row.names=1)
mat <- readRDS("data/geneID_samp_matrix.rds")

# 2. import information about biotype (done on matrix, see 'get_biotypes_bulk.R')

biotype <- read.table(btf, sep=";",header=T)

biotyfreqs <- table(biotype$gene_biotype)

ggplot(tibble("biotype"=biotype$gene_biotype)) + 
  geom_bar(aes(x=forcats::fct_infreq(biotype), y=..count..), fill="darkgreen" ) +
  coord_flip() + labs(title="gene biotypes", x="biotype", y="n genes") + 
  theme_bw()

# for plotting, select only biotypes represented by >10 genes:
tab.bt <- table(biotype$gene_biotype)
ggsel <- biotype$gene_biotype[biotype$gene_biotype %in% names(tab.bt[tab.bt>10])]
mycols <- rep(brewer.pal(8,"Set2"),4)[-c(1,2,3)]
ggbiotypes <- ggplot(tibble("biotype"=ggsel)) + 
  geom_bar( aes(x=forcats::fct_infreq(biotype), y=..count..), fill=mycols) +
  labs(title="gene biotypes (n genes >10)", x="biotype", y="n genes") +
  coord_flip() + theme_bw()
# =========================================================================

# =========================================================================
# Filter matrices by biotype  
# =========================================================================
###### take away all that is not protein-coding from both :
#         -TPM matrix 
#         -raw count matrix
prot.ty <- biotype %>% filter(gene_biotype=="protein_coding")
protcodTPM <- TPM[prot.ty$ensembl_gene_id,]
protcod.mat <- mat[prot.ty$ensembl_gene_id,]
# =========================================================================

# =========================================================================
# observe TPM and raw counts information (only protein coding) 
# =========================================================================
# TPM
# =========================================================================
meanTPM <- apply(protcodTPM,1,mean) # or sumsTPM/dim(TPM)[2]
ggtpm1 <- ggplot(as_tibble(meanTPM),aes(log10(meanTPM+1e-10))) + 
  ggtitle("mean TPM distribution across genes") + 
  geom_histogram(aes(y=..density..), color="white",fill="#dd9c00",alpha=.7) + 
  geom_vline(xintercept=log(1),linetype = "dashed" ) +
  geom_density( color="gold", alpha=.2) +
  theme_bw()

# LETS SIMULATE CANDIDATE THRESHOLDS (logarithmic) and do plot
# Average TPM rule:
logtpmcutoff <- seq(from=-5, to=6, length=100)
cutoff_ngenes <- tibble("logTPM_cutoff" = logtpmcutoff,
                        "how_many_genes" = sapply(logtpmcutoff, function(x){
                          subv <- meanTPM[meanTPM >= 10^x]
                          return(length(subv))  }) #end sapply
)
ggtpm2 <- ggplot(cutoff_ngenes) +
  geom_line(aes(x = logTPM_cutoff, y=how_many_genes)) +
  labs(subtitle = "logTPM threshold vs genes preserved") +
  theme_bw()
# at least one half rule:how many genes have >= half TPMs values over threshold
nolog <- logtpmcutoff^10
cutoff_ngenes$at_least_half <- sapply(nolog, function(x){
  sum(apply(protcodTPM > x, 1 , mean) >= 0.5 )
})

ggtpm3 <-ggplot(cutoff_ngenes) +
  geom_line(aes(x=logTPM_cutoff,y=at_least_half)) +
  labs(subtitle = "at least one half rule") +
  theme_bw()

#save( meanTPM, file=rdatafile)
#load("onlyprot_TPMinf.RData")
###### end TPM information only protein coding
# =========================================================================

# observe RAW COUNTS INFO
# =========================================================================
mean_counts <- apply(protcod.mat, 1, mean)
variance_counts <- apply(protcod.mat, 1, var)

df_counts <- data.frame("mean"=mean_counts, 
                        "variance"=variance_counts)

## fitting data to Neg binomial model: (uses MASS library)
modNB <- glm.nb(mean_counts~variance_counts, data=df_counts)
summary(modNB)
# compare with other poisson
mod <- glm(mean_counts~variance_counts,family="quasipoisson", data=df_counts)
summary(mod)
logLik(mod)# NA
logLik(modNB) # -129917.5 (df=3)
qchisq(0.95,df.residual(modNB))
deviance(modNB)
v = 1/modNB$theta
simx = seq(0.0, 3, 0.1)
gd = data.frame(simx, g = dgamma(simx, shape = 1/v, scale = v))#ggplot(gd,aes(g,simx)) + geom_line()
pred.v = predict(modNB)
pred.d <- cut(predicted,breaks=unique(quantile(predicted,seq(0,100,5)/100)))
m <- tapply(seq(1:21930), pred.d, mean)
v <- tapply(seq(1:21930),pred.d,var)
plot(m,v)
x <- seq(1:21930) # plot(x, x*(1+x/modNB$theta), col="red")
# use the parameter theta to simulate means:
df_counts$simulvar = df_counts$mean*(1+df_counts$mean/modNB$theta)
df_counts$simulvar = df_counts$mean*(1+df_counts$mean/modNB$theta)
ggrawcounts <- ggplot(df_counts, aes(x=log10(mean_counts), y=log10(variance_counts))) +
  geom_point( size=.2, alpha=.2,
              color="darkgreen") + 
  geom_line(aes(log10(mean_counts), log10(simulvar)), color="salmon") + 
  ggtitle("Raw counts (only protein coding)") + theme_bw()
# =========================================================================


# save protcod matrices
saveRDS(protcod.mat, file=out.counts)
saveRDS(protcodTPM, file=out.tpm)

#### print plots to pdf
pdf(out1.pdf)
ggbiotypes
plot_grid(ggtpm1,ggtpm2, ggtpm3, ggrawcounts,ncol=2)
dev.off()
# ======================================================================
# NEW !!! another filtering strategy:
# ======================================================================
# NOTE: manually saved figure as provisoire.pdf
# TODO : fix all figures for this script:
#     - put mean counts by mean variance alone alongside biotypes,
#     - and add line y=2x or fit linear
#     -in second page: only TPM data,
#           add meanTPM by day by age, AND include this last one  : 
EV = 1e-1
LprotcodTPM <- log10(protcodTPM+(1*EV))
keep <- apply(LprotcodTPM, 1, function(row) sum(row > log10(EV)) == length(row))
fi_protcodTPM <- protcodTPM[keep,]
filtered_meanTPM <- apply(fi_protcodTPM,1,mean)
ggplot(as_tibble(filtered_meanTPM),aes(log10(value+(1*EV)))) + 
  ggtitle("Filtered mean TPM distribution across genes") + 
  geom_histogram(aes(y=..density..), color="white",fill="#dd9c00",alpha=.7,
                 bins=60) + 
  geom_density( color="gold", alpha=.2) +
  theme_bw() + labs(x=paste0("Log10(TPM+",EV,")"),
    caption =paste0("Retaining only log10(values) superior to min cutoff: \n",
    "keep <- apply(LogTPM,1, function(x) sum(x > ", EV,") == length(x))"))