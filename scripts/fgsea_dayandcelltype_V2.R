#  Runs GSEA by day and celltype
# from 'static' DEGs
# Saves rds objects into exam_INTER_conditions/static/GSEA/rds/
# * Major modification : only REACTOME pathways tested *
# note: filtered rds file contains only top 15 up and top 15 down pathways
# --
# johaGL
library(fgsea)
library(tidyverse)
library(msigdbr)
library(cowplot)
setwd("~/BulkAnalysis_plusNetwork/")
# take top 1000 ranked genes by pvalue (500 up, 500 down)
# and perform gsea (fgsea)

odir = "exam_INTER_conditions/static/"
# degfile = "shot_dataframe_softfilter.csv"   # this is soft filtered version
# DEdf <- read.table(paste0(odir, degfile), header=T, sep='\t')
DEdf <- readRDS(paste0(odir,"rds/shot_rds_full.rds"))
summary(DEdf)

if (!"symbol" %in% colnames(DEdf)){
  genes_df <- read.table("data/genesinfo.csv",sep="\t",header=T)
  DEdf$symbol <- genes_df[match(DEdf$id, genes_df$Geneid),]$symbol
  head(DEdf)
  print("added symbols to dataframe")
}else{ print("symbol already in dataframe")}

# ======= split by day and keep in a list of dataframes (each df ~1000 genes): =====
DE_l = list()  #  
for (k in c('D0','D2', 'D4', 'D7')){
  tmp <- DEdf %>% filter(day == k)  %>% filter(!is.na(padj)) %>%  ## NEW !isna
    mutate(sens= ifelse(log2FoldChange < 0,"down", "up"))
  # keep top genes by abslfc and pval
  tmp_up <- tmp %>% filter(sens == "up") %>% group_by(type) %>%
    arrange(padj, .by_group = TRUE) %>% slice_min(padj, n=500, with_ties=F) %>%
    filter(abs(log2FoldChange) >= 0.3)
  tmp_down <- tmp %>% filter(sens == "down") %>% group_by(type)   %>%
    arrange(padj,  .by_group = TRUE) %>% slice_min(padj, n=500, with_ties=F) %>%
    filter(abs(log2FoldChange) >= 0.3)
  DE_l[[k]] <-  rbind(tmp_up,tmp_down)
}
infoinput <- data.frame("day_celltype" = c(), "inputsize"=as.integer(c()), 
                        "maxpadj" = c(), "minabslfc" = c())
for (k in c('D0', 'D2', 'D4', 'D7')){
  print("")
  print(paste(":::::",k,":::::"))
  cts <- unique(DE_l[[k]]$type) 
  tmpplots_ <- list()
  for (CT in cts) {
    print(paste("   --->", CT))
    here.df <- DE_l[[k]] %>% filter(type == CT) %>% arrange(desc(log2FoldChange)) %>%
      select(log2FoldChange, symbol, padj)
    nbg = dim(here.df)[1]
    maxp = max(here.df$padj)
    minalf = min(abs(here.df$log2FoldChange))
    infoinput <- rbind(infoinput, c(paste0(k,"_",CT),
                                    nbg, maxp, minalf))
  } }
colnames(infoinput) <- c("day_celltype" , "inputsize", 
                         "maxpadj" , "minabslfc" )
infoinput$inputsize <- as.integer(infoinput$inputsize)
write.table(infoinput, paste0(odir,"GSEA/csv/infoinput.csv"), sep='\t', 
            col.names = T, row.names = F)

plots_ <- list()
for (k in c('D0', 'D2', 'D4', 'D7')){
  cts <- unique(DE_l[[k]]$type) 
  tmpplots_ <- list()
  for (CT in cts) {
    here.df <- DE_l[[k]] %>% filter(type == CT) %>% arrange(desc(log2FoldChange)) %>%
      select(log2FoldChange, symbol)
    gseagenes = here.df %>% pull(log2FoldChange)
    names(gseagenes) = here.df$symbol
    barplot(gseagenes, main=paste(k, CT))
    # same barplot but on ggplot: 
    here.df <- here.df %>% mutate(seq=1:n())
    new <- c(here.df[1,]$symbol) # new labels, first element
    saut = 10
    if (k == 'D2' & (CT == 'ECs' | CT == 'M1')){
      addsum = 10
    }else{addsum = 40}
    for (i in 2:dim(here.df)[1]){
      if( i == saut){
        new = c(new, here.df[i,]$symbol )
        saut = saut + addsum
      }else{ new = c(new, "") } } # end for
    gplo <- ggplot(here.df , aes(x=reorder(symbol, log2FoldChange),y=log2FoldChange)) +
      geom_col() + 
      theme( axis.ticks.x = element_blank(),
             axis.text.x = element_text(angle = 90)) +
      scale_x_discrete(labels= new) + 
      labs(title=paste(k, CT), x = "genes") 
    tmpplots_[[CT]] <- gplo
  }
  plots_[[k]] <- tmpplots_
}
library(cowplot)
pdf(paste0(odir, "GSEA/preGSEA_byday_bycelltype_V2.pdf"), width=12, height=30)
plots_[['D0']][["e1"]] <- NULL
plot_grid(
  plot_grid(plotlist = plots_[['D0']], nrow=2, ncol=3),
  plot_grid(plotlist = plots_[['D2']], nrow= 2, ncol = 3),
  plot_grid(NULL),
  plot_grid(plotlist = plots_[['D4']], nrow=2, ncol = 3),
  plot_grid(NULL),
  plot_grid(plotlist = plots_[['D7']], nrow=2, ncol = 3),
  nrow = 6, rel_heights = c(3,3,1,3,1,3)
)
dev.off()

# ======================== RUNNING GSEA =======================================
thegmt <- msigdbr(species = "Mus musculus", 
                  category = 'C2', 
                  subcategory=c('CP:REACTOME'))
msigdbr_list = split(x = thegmt$gene_symbol, f = thegmt$gs_name)

pathsFull_l <- list("D0"=list(),"D2"=list(),"D4"=list(),"D7"=list())
pathsFiltered  <- list("D0"=list(),"D2"=list(),"D4"=list(),"D7"=list())
for (k in c('D0', 'D2', 'D4', 'D7')){
  cts <- unique(DE_l[[k]]$type) 
  dd <- list()
  dd_f <- list()
  for (CT in cts){
    set.seed(42)
    here.df <- DE_l[[k]] %>% filter(type == CT) %>% arrange(desc(log2FoldChange)) %>%
      select(log2FoldChange, symbol)
    print(c(k, CT))
    gseagenes = here.df %>% pull(log2FoldChange)
    names(gseagenes) <- here.df$symbol
    pathsFull_l[[k]] = list()
    pathsFiltered[[k]] = list()
    print("running fgseaMultilevel, nperm not needed")
    fgseaRes <- fgsea::fgsea(pathways = msigdbr_list, 
                             stats = gseagenes,
                             minSize=3,
                             maxSize=Inf,
                             eps=0) 
    topPathwaysUp <- fgseaRes[ES>0][head(order(padj),n=15), ]
    topPathwaysDown <- fgseaRes[ES<0][head(order(padj), n=15), ]
    combipath <- rbind(topPathwaysUp, topPathwaysDown)
    combipath <- combipath %>% 
      mutate(sens = ifelse(ES > 0, "UP", "DOWN" ))
    dd[[CT]] <- fgseaRes
    dd_f[[CT]] <- combipath
  }
  pathsFull_l[[k]] <- dd
  pathsFiltered[[k]] <- dd_f
}
saveRDS(pathsFull_l, file=paste0(odir,"GSEA/rds/fgseaByDay_fullV2.rds" ))
saveRDS(pathsFiltered, file=paste0(odir, "GSEA/rds/fgseaByDay_filteredV2.rds" ))
print("this Version 2 (V2) yields padj not significant because input was smaller than in first version (no filename suffix)")

