# joha GL 2021
library(dplyr)
library(tidyverse)
library(ggplot2)
library(cowplot)
library(forcats)
library(RColorBrewer)
library(MASS)
library(pheatmap)
library(ggsci) # publishing palettes
library(cowplot)
library(gridExtra)
library(reshape2)
library(ggthemes)
library(scales) # 'viridis_pal' et autres pal
library(ggrepel) # for labels to points
register(MulticoreParam(4)) # TODO:  set n of cores depending of available
library(org.Mm.eg.db)
library(AnnotationDbi)

setwd("~/BulkAnalysis_plusNetwork")
resdir = "plotsDE/"
DEdir = "signaturestypes/"

# https://cran.r-project.org/web/packages/msigdbr/vignettes/msigdbr-intro.html
library("msigdbr") # install.packages("msigdbr")

all_gene_sets = msigdbr(species = "Mus musculus")
setsli <- list()
setsli[["H"]] = msigdbr(species = "Mus musculus", category = "H")
setsli[["C2"]] = msigdbr(species = "Mus musculus", category = "C2")

saveRDS(setsli, file=paste0(DEdir,"DATABASE.rds"))

View(msigdbr_collections())


age = "Old"
vrs = "vA"  # TODO: file version to open
df <- read.table(paste0("signaturestypes/edger_dynINTRA_",age,vrs,".txt"), 
             sep='\t', header=T)


ens2entrez <- function(vectorensemblids){
  res <- AnnotationDbi::select(org.Mm.eg.db,
                               key=vectorensemblids,
                               columns="ENTREZID", 
                               keytype = "ENSEMBL")
  return(res) # a dataframe columns ENSEMBL and ENTREZID
}

haha <- ens2entrez(c("ENSMUSG00000001305", "ENSMUSG00000001143"))

ages=c("Young","Old")
# kegg by pairwise: 
for (ag in ages[1]){
  resu_l <- list()
  gokegg <- list()
  for (ct in allct[3]){
    agectmx <- fmat[,str_detect(colnames(fmat), ag) & str_detect(colnames(fmat), ct)]
    acmeta <- metadata %>% filter(age==ag & type==ct)
    timeps <- sort(unique(acmeta$time))
    Ntp = length(timeps)   #  this variable is used downstream in this loop
    print(timeps)
    timepoints <- factor(acmeta[match(colnames(agectmx),acmeta$newname),]$time,
                         levels=timeps)
    # avoid multiple ensembl to 1 common entrez :
    ens_entrez <- ens2entrez(rownames(agectmx)) %>% as_tibble() %>%
      group_by(ENTREZID) %>% filter(!(n_distinct(ENSEMBL) >= 2))
    # delete rows not being maped to entrezid
    agectmx <- agectmx[ens_entrez$ENSEMBL,]
    rownames(agectmx) <- ens_entrez[match(rownames(agectmx),
                                          ens_entrez$ENSEMBL),]$ENTREZID
    y <- DGEList(counts=agectmx, group=timepoints)
    keep <- filterByExpr(y)
    y <- y[keep, , keep.lib.sizes=FALSE]
    y <- calcNormFactors(y)
    basecontrast <- rep(0,Ntp)
    design <- model.matrix(~0+timepoints)
    y <- estimateDisp(y,design)
    fit <- glmQLFit(y, design)
    tmp_l <- list()
    go_l <- list()
    ke_l <- list()
    for (k in 2:Ntp){
      j <- k-1
      basecontrast[j] <- -1
      basecontrast[k] <- 1
      print(basecontrast)
      qlf <- glmQLFTest(fit, contrast=basecontrast)
      keg <- kegga(qlf,species="Mm")
      go <- goana(qlf, species="Mm")
      gr <- topGO(go)
      kr <- topKEGG(keg)
      dynamicN <- ifelse(Ntp==2, 300 ,ifelse(Ntp==3, 150, 80))
      tr = as.data.frame(topTags(qlf, sort.by="logFC",p.value=0.3, n=dynamicN))
      tr$type = ct
      thiscontrast = paste0(c(timeps[k],timeps[j]),collapse="vs")
      tr$contrast = thiscontrast
      tr$ensemblid = rownames(tr)
      tmp_l[[thiscontrast]] <- tibble(tr)
      go_l[[thiscontrast]] <- gr
      ke_l[[thiscontrast]] <- kr
      basecontrast <- rep(0,Ntp) # reinitialize base contrast
    }
    gokegg[["go"]] <- ke_l
    gokegg[["kegg"]] <- go_l
    resu_l[[ct]] <- bind_rows(tmp_l)
  }#end for ct
  kkk <- bind_rows(resu_l)
  write.table(kkk, file=paste0(DEdir,"edger_dynINTRA_",ag,vrs,".txt"), sep='\t',
              col.names=T, row.names=F)
}#end for
