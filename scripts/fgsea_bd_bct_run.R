#  Runs GSEA by day and celltype
# from 'static' DEGs and non DEGs
# Saves rds objects into exam_INTER_conditions/static/GSEA_bd_bct/rds/
# * Major modification : only REACTOME pathways tested *
# note: filtered rds file contains only top 15 up and top 15 down pathways
# --
# johaGL
library(fgsea)
library(tidyverse)
library(msigdbr)
library(cowplot)
setwd("~/BulkAnalysis_plusNetwork/")
# take top M ranked genes by pvalue (M/2 up, M/2 down) and perform gsea (fgsea)
M = 1000
daysv = c('D0', 'D2', 'D4', 'D7')
DEfile = "rds/shot_rds_full.rds"
odir = "exam_INTER_conditions/static/"
infoinputfile = "GSEA_bd_bct/csv/infoinput"  # csv or md_txt extensions (lines 67 71 )
pdfplotfile = "preGSEA_byday_bycelltype.pdf"
gseards = "GSEA_bd_bct/rds/"
gseaoutfull = "fgsea_bd_bct_full.rds" 
gseaoutfiltered = "fgsea_bd_bct_filtered.rds"

DEdf <- readRDS(paste0(odir, DEfile))
summary(DEdf)
system(paste0("cd ",odir, "; if [ ! -d ", 
              gseards,"   ] ;then mkdir -p ", gseards, ";fi"))

if (!"symbol" %in% colnames(DEdf)){
  genes_df <- read.table("data/genesinfo.csv",sep="\t",header=T)
  DEdf$symbol <- genes_df[match(DEdf$id, genes_df$Geneid),]$symbol
  head(DEdf)
  print("added symbols to dataframe")
}else{ print("symbol already in dataframe")}

# ======= split by day and keep in a list of dataframes (each df M genes): =====
DE_l = list()  #  
for (k in daysv){
  tmp <- DEdf %>% filter(day == k)  %>% 
     mutate(sens= ifelse(log2FoldChange < 0,"down", "up"))
  # keep top genes by abslfc and pval
  tmp_up <- tmp %>% filter(sens == "up") %>% group_by(type) %>%
    arrange(padj, .by_group = TRUE) %>% slice_min(padj, n=M/2, with_ties = F)
  tmp_down <- tmp %>% filter(sens == "down") %>% group_by(type)   %>%
    arrange(padj,  .by_group = TRUE) %>% slice_min(padj, n=M/2, with_ties = F)
  DE_l[[k]] <-  rbind(tmp_up, tmp_down)
}

saveinfoinput <- function(DE_l){
  infoinput <- data.frame("day_celltype" = character(), "inputsize"=numeric(), 
                          "maxpadj" = numeric(), "minabslfc" = numeric())
    for (k in daysv){
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
      infoinput$inputsize <- as.numeric(infoinput$inputsize)
      write.table(infoinput, paste0(odir, infoinputfile, ".csv"), sep='\t', 
                  col.names = T, row.names = F)
      write.table(infoinput %>% mutate(day_celltype = paste("|  ", day_celltype)) %>%
                    mutate(minabslfc = paste(minabslfc, " |")), 
                  file = paste0(odir, infoinputfile, "_md.txt"), sep=" | ", 
                  col.names = T, row.names = F)
    }
saveinfoinput(DE_l)

# ======================== perform Gsea on top ranked genes ====================

# visualize and save barplots before running fgsea
dobarplotbefore <- function(DE_l){
  plots_ <- list()
  for (k in daysv){
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
      saut = 40
      addsum = 40
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
  pdf(paste0(odir, "GSEA_bd_bct/", pdfplotfile), width=12, height=30)
  plots_[['D0']][["e1"]] <- NULL
  print( plot_grid(
    plot_grid(plotlist = plots_[['D0']], nrow=2, ncol=3),
    plot_grid(plotlist = plots_[['D2']], nrow= 2, ncol = 3),
    plot_grid(NULL),
    plot_grid(plotlist = plots_[['D4']], nrow=2, ncol = 3),
    plot_grid(NULL),
    plot_grid(plotlist = plots_[['D7']], nrow=2, ncol = 3),
    nrow = 6, rel_heights = c(3,3,1,3,1,3)
  ) ) # end print
  dev.off()
}
# dobarplotbefore(DE_l)

# ======================== RUNNING GSEA =======================================
thegmt <- msigdbr(species = "Mus musculus", 
                  category = 'C2', 
                  subcategory=c('CP:REACTOME'))
msigdbr_list = split(x = thegmt$gene_symbol, f = thegmt$gs_name)

pathsFull_l <- list("D0"=list(),"D2"=list(),"D4"=list(),"D7"=list())
pathsFiltered  <- list("D0"=list(),"D2"=list(),"D4"=list(),"D7"=list())
for (k in daysv){
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
                             maxSize=Inf) 
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

saveRDS(pathsFull_l, file=paste0(odir,gseards, gseaoutfull ))
saveRDS(pathsFiltered, file=paste0(odir,gseards, gseaoutfiltered  ))


# ================== printing GSEA min padj results ==========================
XX = readRDS(paste0(odir,gseardsl, gseaoutfull))
padjofpaths = array(NA, dim=c(4,6))
rownames(padjofpaths) = daysv
colnames(padjofpaths) = names(XX[['D2']])
for (d in daysv){
  for (n in names(XX[[d]])){
    print(n)
    print(min(XX[[d]][[n]]$padj))
    padjofpaths[d, n] <- min(XX[[d]][[n]]$padj)
  }
}
# ===============================end ========================================
