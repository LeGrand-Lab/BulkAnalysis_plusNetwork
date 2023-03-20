#  Runs GSEA by day and celltype
# from 'static' and 'dynamic' DEGs and non DEGs
# List gene order by padj * sign log2foldchange
# Saves rds objects into exam_INTER_conditions/static/GSEA/rds/
# * Major modification : only REACTOME pathways tested *
# note: filtered rds file contains only top 15 up and top 15 down pathways for shiny
# --
# johaGL + MoullePauline
library(fgsea)
library(ComplexHeatmap)
library(msigdbr)
library(cowplot)
library(treemap)
library(networkD3)
library(dplyr)
library(ggplot2)
library(hrbrthemes)
library(RColorBrewer)
display.brewer.all(colorblindFriendly = T)
library(webshot)
library(htmlwidgets)
library(patchwork)
library(facetious)
library(gridExtra)
library(kableExtra)

##################################
#### Path files to load or created
##################################


setwd("~/BulkAnalysis_plusNetwork/")
#Directories
odir <- "exam_INTER_conditions/"
staticD = paste0(odir,"static/")
gseardsFull = "GSEA/"
for ( sd in c("TableGSEA/","PlotsGSEA/","GSEAanalysis/") ){
      if (dir.exists(paste0(staticD,gseardsFull)) == F) {
        dir.create(paste0(staticD,gseardsFull))
      }
      if (dir.exists(paste0(staticD,gseardsFull,sd)) == F) {
        dir.create(paste0(staticD,gseardsFull,sd))
      }
}

#Loads
DEG_table = "TableDEG/DEG_table_"
#DEfilestatic = paste0(staticD,DEG_table, "[Young_|Old_|][static|dynamic]_full.rds")

#Created
#Table GSEA -> paste0(staticD,TableGSEA, "[Young_|Old_|][static|dynamic]_[full|solftfilter].rds"
#softfilter -> top padj 15 UP + DOWN -> use for heatmaprecap
#Plots -> Heatmap
#Analysis -> report
GSEA_table = paste0(gseardsFull,"TableGSEA/GSEA_table_")
GSEA_plot= "PlotsGSEA/"
  

  #####
  # GSEA analysis
  # static old vs young
  # by cell-type by day
  #####
  DEdf<-readRDS(paste0(staticD,DEG_table,"static_full.rds"))
  days<-unique(DEdf$day)
  typesv<-lapply(days , function(d) unique(DEdf[DEdf$day == d,]$type) )
  names(typesv)<-days
  daysv<-list()
  daysv<-lapply( unique(DEdf$type), function(t) unique(DEdf[DEdf$type == t,]$day) )
  names(daysv)<-unique(DEdf$type)
  Typecellv<-unique(DEdf$type)
  
  
  
  if (!"symbol" %in% colnames(DEdf)){
    genes_df <- read.table("data/genesinfo.csv",sep="\t",header=T)
    DEdf$symbol <- genes_df[match(DEdf$id, genes_df$Geneid),]$symbol
    head(DEdf)
    print("added symbols to dataframe")
  }else{ print("symbol already in dataframe")}
  
  # ======================== RUNNING GSEA =======================================
  thegmt <- msigdbr(species = "Mus musculus", 
                    category = 'C2', 
                    subcategory=c('CP:REACTOME'))
  msigdbr_list = split(x = thegmt$gene_symbol, f = thegmt$gs_name)
  
  List_GSEA_result <- list("D0"=list(),"D2"=list(),"D4"=list(),"D7"=list())
  List_GSEA_result_filt  <- list("D0"=list(),"D2"=list(),"D4"=list(),"D7"=list())
  table_GSEA_result <- data_frame()
  
  for (d in days){
    dd_rank_full <- list()
    dd_rank_filtered <- list()
    for (CT in typesv[[d]]){
      set.seed(42)
      # analyse with signed log2foldchange * -log10(padj)
      here.df_arrange_padj_full<- DEdf %>% filter(type == CT, day == d) %>% drop_na(padj) %>% mutate(rank = log2FoldChange/abs(log2FoldChange) * (-1*log10(padj))) %>% dplyr::select(rank, symbol) %>% arrange(rank)
      print(c(d, CT))
      gseagenes_rank_full = here.df_arrange_padj_full %>% pull(rank)
      names(gseagenes_rank_full) <- here.df_arrange_padj_full$symbol
      List_GSEA_result[[d]] = list()
      List_GSEA_result_filt[[d]] = list()
      print("running fgseaMultilevel, nperm not needed")
      fgseaRes_rank_full <- fgsea::fgsea(pathways = msigdbr_list, 
                                                   stats = gseagenes_rank_full,
                                                   minSize=3,
                                                   maxSize=Inf) 
      topPathways_rank_Up <- fgseaRes_rank_full[ES>0][head(order(padj),n=15), ]
      topPathways_rank_Down <- fgseaRes_rank_full[ES<0][head(order(padj), n=15), ]
      combipath_rank_full <- rbind(topPathways_rank_Up, topPathways_rank_Down)
      combipath_rank_full <- combipath_rank_full %>% 
        mutate(sens = ifelse(ES > 0, "UP", "DOWN" ))
      dd_rank_full[[CT]] <- fgseaRes_rank_full
      dd_rank_filtered[[CT]] <- combipath_rank_full
      fgseaRes_rank_full$day <- d
      fgseaRes_rank_full$type <- CT
      table_GSEA_result<-rbind(table_GSEA_result,fgseaRes_rank_full)
      
    }
    List_GSEA_result[[d]] <- dd_rank_full
    List_GSEA_result_filt[[d]] <- dd_rank_filtered
    
  }
  saveRDS(List_GSEA_result, file=paste0(staticD,GSEA_table, "static_full.rds" ))
  saveRDS(List_GSEA_result_filt, file=paste0(staticD,GSEA_table, "static_softfilter.rds" ))
  write.table(table_GSEA_result %>% mutate(leadingEdge=str_replace_all(as.character(leadingEdge), "[c()\",]" ,"" )), file=paste0(staticD,GSEA_table, "static_full.csv" ), row.names = F)
GSEAsigni<- table_GSEA_result %>% filter(padj<=0.05)

