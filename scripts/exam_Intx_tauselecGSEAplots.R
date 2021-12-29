# Do plots with GSEA_Intx results, the classical ones 

#######
# johaGL
#####

setwd("~/BulkAnalysis_plusNetwork/")
library(tidyverse)
library(DESeq2)
library(cowplot)
library(fgsea)
library(msigdbr)

odir = "exam_INTER_conditions/static/"
resEnsDE <- readRDS(paste0(odir, "rds/Intx_shot_onTauExtract.rds"))
Reac_gmt <- msigdbr(species = "Mus musculus", category = 'C2', subcategory=c('CP:REACTOME'))
msigdbr_list = split(x = Reac_gmt$gene_symbol, f = Reac_gmt$gs_name)
#    "Intx_fgsea_fullNEW.rds"  
fgseaR_ = readRDS(paste0(odir, "GSEA_Intx/rds/", "Intx_fgsea_full.rds" ))

getgeneslist <- function(gseadatafr, d){
  outi <- c()
  for (i in gseadatafr$leadingEdge){
    outi <- c(outi, i)
  }
  moo <- unique(outi)
  tmpdfdeg <- resEnsDE[[d]]
  lfcs <- tmpdfdeg[match(moo, tmpdfdeg$symbol),]$log2FoldChange
  names(lfcs) <- moo
  return(lfcs)
}

if (tablesneeded){
  system(paste0("cd ",odir,"GSEA_Intx/; ", 
                "if [ ! -d csv/ ]; then mkdir csv; else echo 'csv/ exists, no action';fi"))
  print(" ** see keys in full results rds list of lists ** ")
  for (d in names(fgseaR_)){
    print(paste("\n^^^^^", d, "^^^^^"))
    tmpdf <- fgseaR_[[d]]
    tmpdf <-  tmpdf %>% mutate(gene_symbols = sapply(leadingEdge, 
                                                     function(x) paste(x, collapse=", ")),
                               .before = leadingEdge) %>%
      select(-leadingEdge)
    write.table(tmpdf, paste0(odir,"GSEA_Intx/csv/Intx_",d,"_pathways.csv"), sep="\t", 
                col.names = T, row.names=F)
  }
}
plotme <- function(fgseaR_, d){
  fgseaRes <- fgseaR_[[d]]
  summary(fgseaRes$padj)
  NumP = 15
  topPathwaysUp <- fgseaRes[ES>0][head(order(padj),n=NumP), ]
  topPathwaysDown <- fgseaRes[ES<0][head(order(padj), n=NumP), ]
  topPathBoth <- c(topPathwaysUp$pathway, rev(topPathwaysDown$pathway))
  a <- getgeneslist(topPathwaysUp, d)
  b <- getgeneslist(topPathwaysDown, d)
  gseagenes <- c(a, b)
  print(gseagenes)
  ouif = fgsea::plotGseaTable(msigdbr_list[topPathBoth], gseagenes, fgseaRes, 
                              gseaParam = 0.5 , render=F) 
 
  plotsenrichu_ <- list()
  for (i in 1:NumP){
    pup = topPathwaysUp[i,]
    tmpup <- fgsea::plotEnrichment(msigdbr_list[[pup$pathway]], gseagenes) + 
      labs(title= d, subtitle = pup$pathway, 
           caption=paste( "(NES:", pup$NES, ", padj :", pup$padj,")"), render = F)
    if (dim(tmpup$data)[1] <= 2 ){
      plotsenrichu_[[i]] <- NULL
    }else{ plotsenrichu_[[i]] <- tmpup}
  }
  plotsenrichdw_ <- list()
  for (i in 1:NumP){
    pdw = topPathwaysDown[i,]
    tmpdw <- fgsea::plotEnrichment(msigdbr_list[[pdw$pathway]], gseagenes) + 
      labs(title= d, subtitle = pdw$pathway,
           caption=paste( "(NES:", pdw$NES, ", padj :", pdw$padj,")") , render = F)
    if (dim(tmpdw$data)[1] <= 2 ){
      plotsenrichdw_[[i]] <- NULL
    } else {plotsenrichdw_[[i]] <- tmpdw }
  }
  return(list(ouif, plotsenrichu_, plotsenrichdw_))
}


pgrids <- plotme(fgseaR_, "D0") 
pdf(paste0(odir, "GSEA_Intx/", "intx_NEW", "D0",".pdf" ), width=15, height=15)
plot_grid(ggdraw() + draw_label(paste("D0" , ": Top enriched Pathways (GSEA), Old vs Young")),
          plot_grid(NULL, pgrids[[1]] , rel_widths =c(4,7)),
          NA,
          nrow = 3, rel_heights = c(1,9,9))
plot_grid(plotlist = pgrids[[2]], nrow = 5)
plot_grid(plotlist = pgrids[[3]], nrow=5)
dev.off()


pgrids <- plotme(fgseaR_, "D2") 
pdf(paste0(odir, "GSEA_Intx/", "intx_NEW", "D2",".pdf" ), width=15, height=15)
plot_grid(ggdraw() + draw_label(paste("D2" , ": Top enriched Pathways (GSEA), Old vs Young")),
          plot_grid(NULL, pgrids[[1]] , rel_widths =c(4,7)),
          NA,
          nrow = 3, rel_heights = c(1,9,9))
plot_grid(plotlist = pgrids[[2]], nrow = 5)
plot_grid(plotlist = pgrids[[3]], nrow=5)
dev.off()


pgrids <- plotme(fgseaR_, "D4") 
pdf(paste0(odir, "GSEA_Intx/", "intx_NEW", "D4",".pdf" ), width=15, height=15)
plot_grid(ggdraw() + draw_label(paste("D4" , ": Top enriched Pathways (GSEA), Old vs Young")),
          plot_grid(NULL, pgrids[[1]] , rel_widths =c(4,7)),
          NA,
          nrow = 3, rel_heights = c(1,9,9))
plot_grid(plotlist = pgrids[[2]], nrow = 5)
plot_grid(plotlist = pgrids[[3]], nrow=5)
dev.off()


pgrids <- plotme(fgseaR_, "D7") 
pdf(paste0(odir, "GSEA_Intx/", "intx_NEW", "D7",".pdf" ), width=15, height=15)
plot_grid(ggdraw() + draw_label(paste("D7" , ": Top enriched Pathways (GSEA), Old vs Young")),
          plot_grid(NULL, pgrids[[1]] , rel_widths =c(4,7)),
          NA,
          nrow = 3, rel_heights = c(1,9,9))
plot_grid(plotlist = pgrids[[2]], nrow = 5)
plot_grid(plotlist = pgrids[[3]], nrow=5)
dev.off()


### NES vs Tau:
taude <- readRDS(paste0(odir, "rds/TauPlusDEinfo_full.rds"))

nestau <- data.frame("day" = c(), "symbol"=c(),"NES"=c())
for (d in names(fgseaR_)){
  for (r in 1:dim(fgseaR_[[d]])[1]){
    thegenes <- unlist(fgseaR_[[d]][r,]$leadingEdge)
    is.vector(thegenes)
    for (k in 1:length(thegenes)){
      nestau <- rbind(nestau, c(d, thegenes[k],
                                fgseaR_[[d]][r,]$NES))
    }
    
  }
  nestau <- nestau %>% unique()
} 
colnames(nestau) = c("day", "symbol", "NES")

fofo = left_join(nestau, taude, by=c("symbol", "day"))
fofo = fofo %>% select(day, symbol, NES, Tau) %>% filter(!is.na(Tau))
ggplot(fofo, aes(x=NES, y=Tau, color=day)) + 
  geom_point(size=.7, alpha=.6) + scale_color_brewer(palette="Dark2", directionn=-1)




