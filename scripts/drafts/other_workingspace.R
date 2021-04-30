library(dplyr)
library(tidyverse)
library(cowplot)
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
nolog <- logtpmcutoff^10
cutoff_ngenes$at_least_half <- sapply(nolog, function(x){
  sum(apply(TPM > x, 1 , mean) >= 0.5 )
})

ggplot(cutoff_ngenes) +
  geom_line(aes(x=logTPM_cutoff,y=at_least_half)) +
  ggtitle("test") +
  theme_bw()

cutoff_ngenes$otherhalf <- sapply(nolog,function(x){
  sum(apply(TPM,1,function(y) sum(y>=x)))
})


ggtpm3 <- ggplot(cutoff_ngenes) +
  geom_line(aes(x=logTPM_cutoff,y=otherhalf)) +
  ggtitle("manual test") +
  theme_bw()

plot_grid(ggtpm1,ggtpm2, ggtpm3, nrows = 2)

save(cutoff_ngenes, meanTPM,sumsTPM,full_least ,file="bulk_workingspace.RData")

## another version of at_least_half rule:
cutoff_ngenes$at_least_half <- sapply(nolog, function(x){
  sup.or.eq.values <- apply(TPM, 1, function(y) sum(y>=x))
  return(sum(sup.or.eq.values>=127/2))
})
ggplot(cutoff_ngenes) +
  geom_line(aes(x=logTPM_cutoff,y=at_least_half)) +
  ggtitle("at least one half rule") +
  theme_bw()
## ! weird hat !! 


# ====   part pending as cutoff not yet established
cutofftpm <- 
newtpm <- meanTPM[meanTPM>cutofftpm]

ggtpmX <- ggplot(as_tibble(),aes(log10(meanTPM))) + 
  ggtitle(paste0("mean TPM after genes filter, cutoff=",cutofftpm)) + 
  geom_histogram(aes(y=..density..), color="white",fill="#dd9c00",alpha=.7) + 
  geom_density( color="gold", alpha=.2) +
  theme_bw()



# ==================== having fun with simulations
# sim data for plotting stuff
# https://stats.stackexchange.com/questions/32035/checking-poisson-distribution-plot-using-mean-and-variance-relationship
set.seed(17)
mu <- 3
n <- 10
n.trials <- 10000
sim <- replicate(n.trials, { x <- rnbinom(n=n,size=10000,mu=mu); c(mean(x),var(x))})
xy <- apply(sim, 1, function(x) quantile(x, probs=seq(0.001, 0.999, length.out=99)))
plot(xy, xlab="mean", ylab="variance")

for (mu in (1:n-1)){
  sim <- replicate(n.trials, { x <- rnbinom(n=n,size=10000,mu=mu); c(mean(x),var(x))})
  xy <- apply(sim, 1, function(x) quantile(x, probs=seq(0.001, 0.999, length.out=99)))
  lines(xy, col="lightblue")
}


######## ? (from filterTPMandCounts.R)
ggcounts1 <- ggplot(df_counts) +
  geom_point(aes(x=mean_counts, y=variance_counts),size=.2, alpha=.2,
             color="darkgreen") + 
  geom_smooth(data= tibble("mean"=x,"variance"=x*(1+x/modNB$theta)), 
              aes(mean,variance), color="salmon") + 
  scale_y_log10() + scale_x_log10() +
  ggtitle("Raw counts (only protein coding)") + theme_bw()

ggcounts2 <- ggplot(df_counts) +
  geom_point(aes(x=log10(mean_counts), y=log10(variance_counts)),
             size=.2, alpha=.2,
             color="darkgreen") + 
  geom_smooth(data= tibble("mean"=x,"variance"=x*(1+x/modNB$theta)), 
              aes(mean,variance), color="salmon") + 
  ggtitle("Raw counts (only protein coding)") + theme_bw()


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
