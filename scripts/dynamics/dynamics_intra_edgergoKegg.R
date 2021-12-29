# joha GL 2021
library(dplyr)
library(tidyverse)
library(ggplot2)
library(RColorBrewer)
library(MASS)
library(pheatmap)
library(gridExtra)
library(reshape2)
library(ggthemes)
library(scales) # 'viridis_pal' et autres pal
library(ggrepel) # for labels to points
library(org.Mm.eg.db)  
library(AnnotationDbi)
library(edgeR)


setwd("~/BulkAnalysis_plusNetwork")
resdir = "dynintra_edger_extended/"
prefil_cou <- "data/prefiltered_counts.rds"
metadata.rds <- "data/metadata.rds"
# https://cran.r-project.org/web/packages/msigdbr/vignettes/msigdbr-intro.html

vrs = "vB"  # TODO: file version 

ens2entrez <- function(vectorensemblids){
  res <- AnnotationDbi::select(org.Mm.eg.db,
                               key=vectorensemblids,
                               columns="ENTREZID", 
                               keytype = "ENSEMBL")
  return(res) # a dataframe columns ENSEMBL and ENTREZID
}
entrez2symbol <- function(vectorensemblids){
  res <- AnnotationDbi::select(org.Mm.eg.db,
                               key=vectorensemblids,
                               columns="SYMBOL", 
                               keytype = "ENTREZID")
  return(res) # a dataframe columns ENSEMBL and ENTREZID
}
# =============  part Edger with limma kegg and go ===========
# ============================================================ 

fmat <- readRDS(prefil_cou)
metadata <- readRDS(metadata.rds)
metadata <- metadata %>% mutate(timetype=paste0(time,".",type)) 
# rows to keep
keep <- apply(fmat, 1, function(row) ifelse(sum(row >=5)>= 3, TRUE, FALSE) )
fmat <- fmat[keep,]
allct <- sort(unique(metadata$type))
allct <- allct[allct !="Neutro"]  #  "ECs"  "FAPs" "M1"   "M2"   "sCs" 
ages=c("Young","Old")

# kegg by pairwise: 
for (ag in ages){
  resu_l <- list()
  gokegg <- list()
  for (ct in allct){
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
      dynamicN <- ifelse(Ntp==2, 300 ,ifelse(Ntp==3, 150, 80))
      keg <- kegga(qlf,species="Mm")
      go <- goana(qlf, species="Mm")
      gr <- topGO(go, ontology=c("BP","CC", "MF"), 
                   number=dynamicN)
      kr <- topKEGG(keg, number=Inf)
      tr = as.data.frame(topTags(qlf, sort.by="logFC",p.value=0.3, n=dynamicN))
      tr$type = ct
      thiscontrast = paste0(c(timeps[k],timeps[j]),collapse="vs")
      tr$contrast = thiscontrast
      tr$id = rownames(tr)
      tmp_l[[thiscontrast]] <- tibble(tr)
      go_l[[thiscontrast]] <- gr
      ke_l[[thiscontrast]] <- kr
      basecontrast <- rep(0,Ntp) # reinitialize base contrast
    }
    gokegg[[paste0(ct,"_go")]] <- go_l
    gokegg[[paste0(ct,"_kegg")]] <- ke_l
    resu_l[[ct]] <- bind_rows(tmp_l)
  }#end for ct
  kkk <- bind_rows(resu_l)
  write.table(kkk, file=paste0(resdir,"edger_dynINTRA_",ag,vrs,".txt"), sep='\t',
              col.names=T, row.names=F)
  saveRDS(gokegg, file=paste0(resdir,"edger_dynINTRAkeggo",ag,vrs,".rds"))
}#end for
# ============================================================ 
# ============   end part Edger with limma kegg and go


