##
# Connects matrices to biotype csv file  to obtain matrices containing 
# only protein coding genes, and also figures into single pdf: 
#   * `plotsPrelim/pre_filtering_figures.pdf`
# * `data/protcod_counts.rds`
# * `data/protcod_TPM.rds`
##     johaGL 2021

library(dplyr)
library(tidyverse)
library(cowplot)
library(RColorBrewer)
library(MASS)
setwd("~/BulkAnalysis_plusNetwork/")

btf <- "data/biotype_bulk.csv"

#outputs
out.counts <- "data/protcod_counts.rds"
out.tpm <- "data/protcod_TPM.rds"
out1.pdf <- "plotsPrelim/biotypes.pdf"
out1.png <-  "plotsPrelim/biotypes.png"
out2.pdf <- "plotsPrelim/rawcounts_obs.pdf"
out3.pdf <- "plotsPrelim/TPM_obs.pdf"
# =========================================================================
# import matrices and visualize properties pre-filtering
# =========================================================================
# 1. import TPM matrix and raw count matrix:
TPM <- read.table("data/TPMmatrix.csv", header=T,sep="\t",row.names=1)
mat <- readRDS("data/geneID_samp_matrix.rds")

# 2. import information about biotype (done on matrix, see 'get_biotypes_bulk.R')

biotype <- read.table(btf, sep=";",header=T)

biotyfreqs <- table(biotype$gene_biotype)

# for plotting, select only biotypes represented by >10 genes:
tab.bt <- table(biotype$gene_biotype)
data2plot = table(biotype$gene_biotype) %>% data.frame() 
sumFreq=sum(as.numeric(data2plot$Freq))
data2plot$gene_pourcent = round(as.numeric(data2plot$Freq)/sumFreq*100,2)
data2plot = data2plot[-as.numeric(rownames(data2plot[data2plot$Var1 == "NOTFOUND",])), ]

ggbiotypes<-ggplot(dplyr::filter(data2plot,Freq>100), aes(x=reorder(Var1,-gene_pourcent),y=gene_pourcent,fill=Var1)) +
  geom_bar(stat="identity", show.legend = FALSE) + 
  labs(title="Gene biotypes (n genes >100)", x="Biotype", y="Pourcent Genes")+
  coord_flip() + theme_bw()+
  theme (axis.text.y = element_text(color="black",  face="bold",size =12),
         axis.text.x = element_text(color="black",  face="bold",size =10, angle =45))
save_plot(out1.png,ggbiotypes, base_height=8, base_width=5)
# =========================================================================

# =========================================================================
# Filter matrices by biotype  
# =========================================================================
###### take away all that is not protein-coding from both :
#         -TPM matrix 
#         -raw count matrix
prot.ty <- biotype %>% dplyr::filter(gene_biotype=="protein_coding")
protcodTPM <- TPM[prot.ty$ensembl_gene_id,]
protcod.mat <- mat[prot.ty$ensembl_gene_id,]
# =========================================================================
# dim(TPM)
# [1] 55487   127
# > dim(protcodTPM)
# [1] 21930   127
print(paste("we got rid of", (dim(TPM)[1]-dim(protcodTPM)[1]), "genes"))


# EXPLORATIONs
# =========================================================================
# observe  and raw counts and TPM (only protein coding) 
# =========================================================================

#  * RAW COUNTS *
# =========================================================================
mean_counts <- apply(protcod.mat, 1, mean)
variance_counts <- apply(protcod.mat, 1, var)

df_counts <- data.frame("mean"=mean_counts, 
                        "variance"=variance_counts)

# do simple linear model on logarithms (our data)
df_counts <- df_counts %>% mutate(l_m = log10(mean+1),
                                  l_v = log10(variance+1))
simplelinear <- lm( df_counts$l_v ~ df_counts$l_m)
(simplelinear$coefficients)  #  0.1582645     1.8424868 

# note: geom_smooth with "lm" replaces abline(lm()):
ggrawcounts <- ggplot(df_counts, aes(x=log10(mean+1), y=log10(variance+1))) +
  geom_point( size=.2, alpha=.2,
              color="darkgreen") + 
  geom_smooth(method="loess", color="lightblue", se=F, size=.5)   +
  geom_smooth(method="lm",se=T, size=.5, color="salmon", 
              linetype="dashed") +
  labs(title="Raw counts (only protein coding)",
          caption = "blue color : loess regression \n
      salmon dashed : linear regression 'lm(log10var~log10Mean)'") + theme_bw()

# end raw counts
# =========================================================================

# * TPM *
# =========================================================================
# meanTPM gene by gene : 
meanTPM <- apply(protcodTPM,1,mean) # or sumsTPM/dim(TPM)[2]
ggtpm1 <- ggplot(as_tibble(meanTPM),aes(log10(meanTPM+1e-6))) + 
  ggtitle("mean TPM distribution across genes") + 
  geom_histogram(aes(y=..density..), color="white",fill="#dd9c00",alpha=.7) + 
  geom_vline(xintercept=log10(1),linetype = "dashed" , color="darkgray") +
  geom_vline(xintercept=log10(0.1),linetype = "dashed" ,color="darkgray") +
  geom_vline(xintercept=log10(0.01),linetype = "dashed", color="darkgray" ) +
  geom_density( color="gold", alpha=.2) +
  theme_bw()

# LETS SIMULATE CANDIDATE THRESHOLDS (logarithmic) and do plot
# Average TPM rule:
logtpmcutoff <- seq(from=-5, to=6, length=100)
cutoff_ngenes <- tibble("log10TPM_cutoff" = logtpmcutoff,
                        "how_many_genes" = sapply(logtpmcutoff, function(x){
                          subv <- meanTPM[meanTPM >= 10^x]
                          return(length(subv))  }) #end sapply
)
ggtpm2 <- ggplot(cutoff_ngenes) +
  geom_line(aes(x = log10TPM_cutoff, y=how_many_genes)) +
  labs(subtitle = "log10TPM threshold vs genes preserved") +
  theme_bw()
# at least one half rule:how many genes have >= half TPMs values over threshold
nolog <- logtpmcutoff^10
cutoff_ngenes$at_least_half <- sapply(nolog, function(x){
  sum(apply(protcodTPM > x, 1 , mean) >= 0.5 )
})

ggtpm3 <-ggplot(cutoff_ngenes) +
  geom_line(aes(x=log10TPM_cutoff,y=at_least_half)) +
  labs(subtitle = "at least one half rule",
       caption= "*utility is not clear*") +
  theme_bw()

# ======================================================================
# NEW !!! another filtering strategy:
# ======================================================================
# FIXED : fix all figures for this script:
#     -: only TPM data,
#    add meanTPM by day by age, AND include this last one  : 
EV = 1e-1
LprotcodTPM <- log10(protcodTPM+(1*EV))
keep <- apply(LprotcodTPM, 1, function(row) sum(row > log10(EV)) == length(row))
fi_protcodTPM <- protcodTPM[keep,]
filtered_meanTPM <- apply(fi_protcodTPM,1,mean)

ggTPMok <- ggplot(as_tibble(filtered_meanTPM),aes(log10(value+(1*EV)))) + 
  ggtitle("After Filter: mean TPM distribution across genes") + 
  geom_histogram(aes(y=..density..), color="white",fill="#dd9c00",alpha=.7,
                 bins=60) + 
  geom_density( color="gold", alpha=.2) +
  theme_bw() + labs(x=paste0("Log10(TPM+",EV,")"),
           caption =paste0("Retaining only log10(values) superior to min cutoff: \n",
        "keep <- apply(LogTPM,1, function(x) sum(x > ", log10(EV),") == length(x))"))

#  end TPM information only protein coding
# =========================================================================

# save protcod matrices, they are protein BUT NOT FILTERED
# filter ideas are applied later on
# =========================================================================
saveRDS(protcod.mat, file=out.counts)
saveRDS(protcodTPM, file=out.tpm)

#### print plots to pdf
pdf(out1.pdf)
ggbiotypes
dev.off()
pdf(out2.pdf)
plot_grid(ggrawcounts)
dev.off()
pdf(out3.pdf, width= 7)
plot_grid(
  ggTPMok,
  plot_grid(ggtpm1,ggtpm2, ncol=2, labels=c("B","C")),
  nrow=2, labels=c('A'))
dev.off()

###########################################
#END
###########################################
# interesting idea but did not work
## fitting data to Neg binomial model: (uses MASS library)
# https://stats.idre.ucla.edu/r/dae/negative-binomial-regression/
# subdf <- df_counts %>% sample_n(10000) %>% 
#   mutate(meancol=round(mean,0)+1)
# modNB <- glm.nb(meancol ~ 1, data=subdf)
# summary(modNB)
# # compare with other poisson
# mod <- glm(variance~mean ,family="quasipoisson", data=df_counts)
# summary(mod)
# logLik(mod)# NA
# logLik(modNB) # 
# qchisq(0.95,df.residual(modNB))
# deviance(modNB)
# v = 1/modNB$theta
# simx = seq(0.0, 3, 0.1)
# gd = data.frame(simx, g = dgamma(simx, shape = 1/v, scale = v))#ggplot(gd,aes(g,simx)) + geom_line()
# pred.v = predict(modNB)
# pred.d <- cut(pred.v,breaks=unique(quantile(pred.v,seq(0,100,5)/100)))
# msim <- tapply(seq(1:length(pred.d)), pred.d, mean)
# vsim <- tapply(seq(1:length(pred.d)), pred.d, var)
# simregr <- lm(log10(vsim+1)~log10(msim+1))
# plot(msim,vsim )
# abline(simregr)
# fitting did not work ! 
#x <- seq(1:21930); plot(x, x*(1+x/modNB$theta), col="red")
# ggcountssimul <- ggplot(as_tibble(log10(msim+1),log10(vsim+1))) +
#   geom_point(aes(log10(msim+1),log10(vsim+1))) + 
#   geom_abline(slope=simregr$coefficients[[2]],
#               intercept=simregr$coefficients[[1]], color="lightblue") +
#   labs(main="",
#      subtitle = "Simulated data from NegBin distribution",
#         caption=paste(  "Theta :",modNB$theta, 
#                       "obtained from fitting model on realdata"))