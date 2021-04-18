# Joha GL
##
#library(reticulate)
#use_condaenv(condaenv="~/prog_bio/anaconda3/envs/MYVENV",
#             conda="~/prog_bio/anaconda3/bin/conda")
library(dplyr)
library(tidyverse)
library("BiocParallel")
register(MulticoreParam(4)) # TODO:  set n of cores depending of available
library(edgeR)
# INSIDE each age separately, i.e "INTRA" dynamics
# compare dynamics of expression by pairwise between adjacent days
setwd("~/BulkAnalysis_plusNetwork/")
prefil_cou <- "data/prefiltered_counts.rds"
metadata.rds <- "data/metadata.rds"
genes_df <- read.table("data/genesinfo.csv",sep="\t",header=T)

vrs <- "vA"  ## TODO: give name to file version (A,B,C....)

resdir = "plotsDE/"
DEdir = "signaturestypes/"
fmat <- readRDS(prefil_cou)
metadata <- readRDS(metadata.rds)
metadata <- metadata %>% mutate(timetype=paste0(time,".",type)) 

# rows to keep
keep <- apply(fmat, 1, function(row) ifelse(sum(row >=5)>= 3, TRUE, FALSE) )
fmat <- fmat[keep,]

allct <- sort(unique(metadata$type))
allct <- allct[allct !="Neutro"]  #  "ECs"  "FAPs" "M1"   "M2"   "sCs" 

ages=c("Young","Old")
for (ag in ages){
  resu_l <- list()
  for (ct in allct){
    agectmx <- fmat[,str_detect(colnames(fmat), ag) & str_detect(colnames(fmat), ct)]
    acmeta <- metadata %>% filter(age==ag & type==ct)
    timeps <- sort(unique(acmeta$time))
    Ntp = length(timeps)   #  this variable is used downstream in this loop
    print(timeps)
    timepoints <- factor(acmeta[match(colnames(agectmx),acmeta$newname),]$time,
                         levels=timeps)
    y <- DGEList(counts=agectmx, group=timepoints)
    keep <- filterByExpr(y)
    y <- y[keep, , keep.lib.sizes=FALSE]
    y <- calcNormFactors(y)
    basecontrast <- rep(0,Ntp)
    design <- model.matrix(~0+timepoints)
    y <- estimateDisp(y,design)
    fit <- glmQLFit(y, design)
    tmp_l <- list()
    for (k in 2:Ntp){
      j <- k-1
      basecontrast[j] <- -1
      basecontrast[k] <- 1
      print(basecontrast)
      qlf <- glmQLFTest(fit, contrast=basecontrast)
      dynamicN <- ifelse(Ntp==2, 300 ,ifelse(Ntp==3, 150, 80))
      tr = as.data.frame(topTags(qlf, sort.by="logFC",p.value=0.3, n=dynamicN))
      tr$type = ct
      tr$contrast = paste0(c(timeps[k],timeps[j]),collapse="vs")
      tr$ensemblid = rownames(tr)
      tmp_l[[k]] <- tibble(tr)
      basecontrast <- rep(0,Ntp)
    }
    resu_l[[ct]] <- bind_rows(tmp_l)
  }#end for ct
  kkk <- bind_rows(resu_l)
  write.table(kkk, file=paste0(DEdir,"edger_dynINTRA_",ag,vrs,".txt"), sep='\t',
              col.names=T, row.names=F)
}#end for



# EXAMPLE M2: design
# timepointsD2 timepointsD4 timepointsD7
# 1            0            0            1
# 2            1            0            0
# 3            0            1            0
# 4            1            0            0
# 5            0            1            0
# 6            0            0            1
# 7            0            1            0
# 8            1            0            0
# 9            0            0            1
# attr(,"assign")
# [1] 1 1 1
# attr(,"contrasts")
# attr(,"contrasts")$timepoints
# [1] "contr.treatment"