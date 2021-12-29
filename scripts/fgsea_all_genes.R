#  Runs GSEA by day and celltype
# from 'static' DEGs and non DEGs
# Saves rds objects into exam_INTER_conditions/static/GSEA/rds/
# * Major modification : only REACTOME pathways tested *
# note: filtered rds file contains only top 15 up and top 15 down pathways
# --
# johaGL
library(fgsea)
library(tidyverse)
library(msigdbr)
library(cowplot)
library(tidyverse)
library(ComplexHeatmap)
library(msigdbr)
library(cowplot)
library(treemap)
library(networkD3)
library(dplyr)
library(ggplot2)
library(hrbrthemes)
library(webshot)
library(htmlwidgets)
library(gridExtra)
library(kableExtra)
setwd("~/BulkAnalysis_plusNetwork2/")
daysv = c('D0', 'D2', 'D4', 'D7')
DEfile = "shot_rds_full.rds"
odir = "exam_INTER_conditions/static/"
pdfplotfile = "preGSEA_byday_bycelltype.pdf"
gseardsFull = "GSEA/rds/"
gseaoutfull2log2Foldchange = "fgsea_log2FoldChange_full.rds" 
gseaoutfull2padj = "fgsea_padj_full.rds" 
gseaout_filtered_log2Foldchange = "fgsea_log2FoldChange_filtered.rds" 
gseaout_filtered_full2padj = "fgsea_padj_filtered.rds" 

analyse<-"analysis/AnalysisGSEA"

DEdf <- readRDS(paste0(odir, DEfile))
summary(DEdf)

if (!"symbol" %in% colnames(DEdf)){
  genes_df <- read.table("data/genesinfo.csv",sep="\t",header=T)
  DEdf$symbol <- genes_df[match(DEdf$id, genes_df$Geneid),]$symbol
  head(DEdf)
  print("added symbols to dataframe")
}else{ print("symbol already in dataframe")}

# ======= split by day and keep in a list of dataframes (each df M genes): =====
DE_l_full = list()  #  
for (k in daysv){
  tmp <- DEdf %>% filter(day == k)  %>% 
     mutate(sens= ifelse(log2FoldChange < 0,"down", "up"))
  # keep top genes by abslfc and pval
  tmp_up <- tmp %>% filter(sens == "up") %>% group_by(type) %>%
    arrange(padj, .by_group = TRUE)
  tmp_down <- tmp %>% filter(sens == "down") %>% group_by(type)   %>%
    arrange(padj,  .by_group = TRUE)
  DE_l_full[[k]] <-  rbind(tmp_up, tmp_down)
}
# ======================== RUNNING GSEA =======================================
thegmt <- msigdbr(species = "Mus musculus", 
                  category = 'C2', 
                  subcategory=c('CP:REACTOME'))
msigdbr_list = split(x = thegmt$gene_symbol, f = thegmt$gs_name)

pathsFull_log2foldchange_full <- list("D0"=list(),"D2"=list(),"D4"=list(),"D7"=list())
pathsFiltered_log2foldchange_full  <- list("D0"=list(),"D2"=list(),"D4"=list(),"D7"=list())

pathsFull_rank_full <- list("D0"=list(),"D2"=list(),"D4"=list(),"D7"=list())
pathsFiltered_rank_full  <- list("D0"=list(),"D2"=list(),"D4"=list(),"D7"=list())

for (k in daysv){
  cts <- unique(DE_l_full[[k]]$type) 
  dd_log2foldchange_full <- list()
  dd_log2foldchange_filtered <- list()
  
  dd_rank_full <- list()
  dd_rank_filtered <- list()
  
  for (CT in cts){
    set.seed(42)
    # analyse with log2foldchange
    here.df_arrange_log2foldchange_full <- DE_l_full[[k]] %>% filter(type == CT) %>% arrange(desc(log2FoldChange)) %>%
      dplyr::select(log2FoldChange, symbol)
    gseagenes_log2foldchange_full = here.df_arrange_log2foldchange_full %>% pull(log2FoldChange)
    names(gseagenes_log2foldchange_full) <- here.df_arrange_log2foldchange_full$symbol
    pathsFull_log2foldchange_full[[k]] = list()
    pathsFiltered_log2foldchange_full[[k]] = list()
    print("running fgseaMultilevel, nperm not needed")
    fgseaRes_log2foldchange_full <- fgsea::fgsea(pathways = msigdbr_list, 
                             stats = gseagenes_log2foldchange_full,
                             minSize=3,
                             maxSize=Inf) 
    topPathways_log2foldchange_Up <- fgseaRes_log2foldchange_full[ES>0][head(order(padj),n=15), ]
    topPathways_log2foldchange_Down <- fgseaRes_log2foldchange_full[ES<0][head(order(padj), n=15), ]
    combipath_log2foldchange_full <- rbind(topPathways_log2foldchange_Up, topPathways_log2foldchange_Down)
    combipath_log2foldchange_full <- combipath_log2foldchange_full %>% 
      mutate(sens = ifelse(ES > 0, "UP", "DOWN" ))
    dd_log2foldchange_full[[CT]] <- fgseaRes_log2foldchange_full
    dd_log2foldchange_filtered[[CT]] <- combipath_log2foldchange_full
    
    # analyse with signed log2foldchange * -log10(padj)
    here.df_arrange_padj_full<- DE_l_full[[k]] %>% filter(type == CT) %>% drop_na(padj) %>% mutate(rank = log2FoldChange/abs(log2FoldChange) * (-1*log10(padj))) %>% dplyr::select(rank, symbol)
    print(c(k, CT))
    gseagenes_rank_full = here.df_arrange_padj_full %>% pull(rank)
    names(gseagenes_rank_full) <- here.df_arrange_padj_full$symbol
    pathsFull_rank_full[[k]] = list()
    pathsFiltered_rank_full[[k]] = list()
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
    
  }
  pathsFull_log2foldchange_full[[k]] <- dd_log2foldchange_full
  pathsFiltered_log2foldchange_full[[k]] <- dd_log2foldchange_filtered
  pathsFull_rank_full[[k]] <- dd_rank_full
  pathsFiltered_rank_full[[k]] <- dd_rank_filtered
  
}

saveRDS(pathsFull_log2foldchange_full, file=paste0(odir,gseardsFull, gseaoutfull2log2Foldchange ))
saveRDS(pathsFiltered_log2foldchange_full, file=paste0(odir,gseardsFull, gseaout_filtered_log2Foldchange  ))

saveRDS(pathsFull_rank_full, file=paste0(odir,gseardsFull, gseaoutfull2padj ))
saveRDS(pathsFiltered_rank_full, file=paste0(odir,gseardsFull, gseaout_filtered_full2padj  ))

# ================== printing GSEA min padj results ==========================
XX_full = readRDS(paste0(odir,gseardsFull, gseaoutfull2padj))
padjofpathsNumbersigni_full = array(NA, dim=c(4,6))
padjofpathsMin_full = array(NA, dim=c(4,6))
rownames(padjofpathsNumbersigni_full) = daysv
colnames(padjofpathsNumbersigni_full) = names(XX_full[['D2']])
rownames(padjofpathsMin_full) = daysv
colnames(padjofpathsMin_full) = names(XX_full[['D2']])
GSEAconcat<-tibble()
for (d in daysv){
  for (n in names(XX_full[[d]])){
    GSEAconcat<-rbind(GSEAconcat, XX_full[[d]][[n]] %>% mutate(day=d) %>% mutate(type=n))
    print(length(unique(XX_full[[d]][[n]][XX_full[[d]][[n]]$padj <= 0.05,]$pathway)))
    padjofpathsNumbersigni_full[d, n] <-length(unique(XX_full[[d]][[n]][XX_full[[d]][[n]]$padj <= 0.05,]$pathway))
    padjofpathsMin_full[d, n] <-min(XX_full[[d]][[n]]$padj)
  }
}


vectTypeCell <- c("Neutro","M2","M1","ECs","FAPs","sCs")
GSEAsigni<- GSEAconcat %>% filter(padj<=0.05)

PathwayCombiDay<-function(n){
  #keep table with signi gene
  tempoSigniDEpathway<- GSEAsigni %>% filter(type==n)
  #keep only uniqueID sort by days
  uniqueIdByDay = lapply(unique(tempoSigniDEpathway$day), function(x) filter(tempoSigniDEpathway, day == x) %>% dplyr::select(pathway) %>% unlist() %>% as.character())
  names(uniqueIdByDay) = unique(tempoSigniDEpathway$day)
  #found all combinaison where one same gene is on 1 or 2 or 3 or 4 days 
  allcombinaison =ComplexHeatmap::make_comb_mat(uniqueIdByDay)
  return(allcombinaison)
}
PathwayCombiType<-function(n){
  #keep table with signi ene
  tempoSigniDEpathway<- GSEAsigni %>% filter(day==n)
  #keep only uniqueID sort by days
  uniqueIdByDay = lapply(unique(tempoSigniDEpathway$type), function(x) filter(tempoSigniDEpathway, type == x) %>% dplyr::select(pathway) %>% unlist() %>% as.character())
  names(uniqueIdByDay) = unique(tempoSigniDEpathway$type)
  #found all combinaison where one same gene is on 1 or 2 or 3 or 4 days 
  allcombinaison =ComplexHeatmap::make_comb_mat(uniqueIdByDay)
  return(allcombinaison)
}

ScsPathwaycombi<-PathwayCombiDay("sCs")
ECsPathwaycombi<-PathwayCombiDay("ECs")
FAPsPathwaycombi<-PathwayCombiDay("FAPs")
M2Pathwaycombi<-PathwayCombiDay("M2")
M1Pathwaycombi<-PathwayCombiDay("M1")
D0Pathwaycombi<-PathwayCombiType("D0")
D2Pathwaycombi<-PathwayCombiType("D2")
D4Pathwaycombi<-PathwayCombiType("D4")
D7Pathwaycombi<-PathwayCombiType("D7")

vectDay1<-c('D0', 'D2', 'D4', 'D7')
vectDay2<- c( 'D2', 'D4', 'D7')
sourceNode<-c()
targetNode<-c()
valueEdge<-c()
Nbconnection<-c()
for ( typeCell in vectTypeCell ) {
  tempoSigniDEgene<- GSEAsigni %>% filter(type==typeCell)
  uniqueIdByDay = lapply(unique(tempoSigniDEgene$day), function(x) filter(tempoSigniDEgene, day == x) %>% dplyr::select(pathway) %>% unlist() %>% as.character())
  names(uniqueIdByDay) = unique(tempoSigniDEgene$day)
  #UpSetR::upset(UpSetR::fromList(uniqueIdByDay))
  allcombinaison =ComplexHeatmap::make_comb_mat(uniqueIdByDay)
  Sup2Combi<-allcombinaison[comb_degree(allcombinaison) >= 2]
  Node_Target<-names(comb_size(Sup2Combi))
  num_combi<-1
  for ( name in Node_Target){
    particulName<-str_extract_all(name,boundary("character"))
    Nbconnection<-c(Nbconnection,rep(str_count(name,"1"),str_count(name,"1")-1))
    dayDE<-c()
    if ( length(particulName[[1]]) <4 ){ vectDay = vectDay2} else { vectDay = vectDay1 }
    for (i in 1:(length(particulName[[1]])-1)){
      if ( particulName[[1]][[i]] == "1" ){
        dayDE<-c(dayDE,vectDay[i])
        if ( particulName[[1]][[i+1]] == "1" ){
          sourceNode<-c(sourceNode,paste0(vectDay[i],"_",typeCell))
          targetNode<-c(targetNode,paste0(vectDay[i+1],"_",typeCell))
          valueEdge<-c(valueEdge,unname(comb_size(Sup2Combi)[num_combi]))
          dayDE<-c(dayDE,vectDay[i+1])
        }
        else {
          if (i+2 <= length(particulName[[1]]) && particulName[[1]][[i+2]] == "1"){
            sourceNode<-c(sourceNode,paste0(vectDay[i],"_",typeCell))
            targetNode<-c(targetNode,paste0(vectDay[i+2],"_",typeCell))
            valueEdge<-c(valueEdge,unname(comb_size(Sup2Combi)[num_combi]))
            dayDE<-c(dayDE,vectDay[i+2])
          }
          else {
            if(i+3 <= length(particulName[[1]]) && particulName[[1]][[i+3]] == "1"){
              sourceNode<-c(sourceNode,paste0(vectDay[i],"_",typeCell))
              targetNode<-c(targetNode,paste0(vectDay[i+3],"_",typeCell))
              valueEdge<-c(valueEdge,unname(comb_size(Sup2Combi)[num_combi]))
              dayDE<-c(dayDE,vectDay[i+3])
            }
          }
        }
      }
    }
    num_combi=num_combi+ 1
  }
} 

# Make a connection data frame
sourceNode<-c(sourceNode,"D0_M2","D0_M1")
targetNode<-c(targetNode,"D2_M2","D2_M1")
valueEdge<-c(valueEdge,0,0)
Nbconnection<-c(Nbconnection,0,0)
links <- data.frame(
  source=sourceNode,
  target=targetNode,
  value=valueEdge
)
# From these flows we need to create a node data frame: it lists every entities involved in the flow
nodes <- data.frame(
  name=c(as.character(links$source), as.character(links$target)) %>% 
    unique()
)
nodes$group <- as.factor(c(as.character(links$source), as.character(links$target)) %>% 
                           unique())
# With networkD3, connection must be provided using id, not using real name like in the links dataframe.. So we need to reformat it.
links$IDsource <- match(links$source, nodes$name)-1 
links$IDtarget <- match(links$target, nodes$name)-1
links$group <- as.factor(Nbconnection)

# prepare color scale: I give one specific color for each node.
my_color <- 'd3.scaleOrdinal() .domain(["2","3","4","0","D0_ECs", "D2_ECs","D4_ECs","D7_ECs", "D0_FAPs", "D2_FAPs", "D4_FAPs","D7_FAPs","D0_sCs", "D2_sCs","D4_sCs","D7_sCs","D0_M1", "D2_M1","D4_M1","D0_M2" ,"D2_M2", "D4_M2","D7_M2"]) .range(["#D3DDDC","#C7B19C","#446455","white","#0072B2", "#0072B2" , "#0072B2", "#0072B2", "#F0E442", "#F0E442", "#F0E442", "#F0E442","#56B4E9","#56B4E9","#56B4E9","#56B4E9","#D55E00","#D55E00","#D55E00", "#CC79A7", "#CC79A7", "#CC79A7","#CC79A7"])'

# Make the Network. I call my colour scale with the colourScale argument
p <- sankeyNetwork(Links = links, Nodes = nodes, Source = "IDsource", Target = "IDtarget", 
                   Value = "value", NodeID = "name", colourScale=my_color, LinkGroup="group", NodeGroup="group", sinksRight = FALSE,
                   nodePadding=50, nodeWidth = 5,fontSize=10,iterations = 55,margin = list("right"=250), height = 650, width = 850)
p



save_plot(paste0(odir,analyse,"PlotDEGOnlyOneD.png"),PlotDEGonlyOneDay,base_width = 4, base_height = 7)
