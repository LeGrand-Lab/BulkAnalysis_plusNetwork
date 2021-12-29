########################################
#Found genes DE on one, or more days
#Observe the dynamic of expression for Old in comparison to Young across days on one Type cell -> sankeyNetwork
#Export Table with gene ID and a code express Dynamic D0upD2down
#Repeat this operation  
######################################

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
odir <- "exam_INTER_conditions/static/"
resdir <-'analysis/DEGOnOneorMoreDays/'

analyse<-"analysis/DynamiqueDEGbytypeCell/"
pathcombinaisontable<-paste0(odir,analyse,"ListOfDynamicDEofGenes.xlsx")
pathcombinaisontableUPDOWN<-paste0(odir,analyse,"ListOfDynamicUpDownDEofGenes.xlsx")

fullDEsta = readRDS(paste0(odir, "rds/shot_rds_full.rds"))


cellcolors = list("ECs"="#0072B2","FAPs"="#F0E442","M1" = "#D55E00","M2" =  "#CC79A7","Neutro" =  "#009E73","sCs" = "#56B4E9" )

signiDEgene <- fullDEsta %>% filter(padj<=0.05)
vectTypeCell <- c("Neutro","M2","M1","ECs","FAPs","sCs")
vectDay1<-c("D0","D2","D4","D7")
vectDay2<-c("D2","D4","D7")

###  CALCULATE
# Distribution general of genes DE log2foldchange

dataFrameTotal=data.frame(signiDEgene,nbDayDE=rep(0,length(signiDEgene$baseMean)),Dynamic=rep(0,length(signiDEgene$baseMean)))
daystable<-c()
celltypetable<-c()
weightNode<-c()
ViolinDay<-c()
ViolinCell<-c()
ViolinName<-c()
ViolinValue<-c()
for ( typeCell in vectTypeCell ) {
  #keep table with signi gene
  tempoSigniDEgene<- signiDEgene %>% filter(type==typeCell)
  #keep only uniqueID sort by days
  uniqueIdByDay = lapply(unique(tempoSigniDEgene$day), function(x) filter(tempoSigniDEgene, day == x) %>% dplyr::select(id) %>% unlist() %>% as.character())
  names(uniqueIdByDay) = unique(tempoSigniDEgene$day)
  #found all combinaison where one same gene is on 1 or 2 or 3 or 4 days 
  allcombinaison =ComplexHeatmap::make_comb_mat(uniqueIdByDay)
  #keep only gene DE one only one DAY
  Sup2Combi<-allcombinaison[comb_degree(allcombinaison) >= 1]
  
  #prepare vector to plot
  Node_Target<-names(comb_size(Sup2Combi))
  num_combi<-1
  for ( name in Node_Target){
    #Make table
    particulName<-str_extract_all(name,boundary("character"))
    if ( length(particulName[[1]]) <4 ){ vectDay = vectDay2} else { vectDay = vectDay1 }
    for (i in 1:(length(particulName[[1]]))){
      if ( particulName[[1]][[i]] == "1" ){
        dayDE<-vectDay[i]
        daystable<-c(daystable,vectDay[i])
        celltypetable<-c(celltypetable, typeCell)
        weightNode<-c(weightNode,unname(comb_size(Sup2Combi)[num_combi]))
      }
    }
    for (idG in extract_comb(Sup2Combi,name) ){
      dataFrameTotal[dataFrameTotal$id == idG & dataFrameTotal$type == typeCell ,]$nbDayDE<-1
      dataFrameTotal[dataFrameTotal$id == idG  & dataFrameTotal$type == typeCell ,]$Dynamic<-dayDE
      ViolinDay<-c(ViolinDay,dayDE)
      ViolinCell<-c(ViolinCell,typeCell)
      ViolinName<-c(ViolinName,unname(comb_size(Sup2Combi)[num_combi]))
      ViolinValue<-c(ViolinValue,dataFrameTotal[dataFrameTotal$id == idG & dataFrameTotal$type == typeCell & dataFrameTotal$day == dayDE ,]$log2FoldChange)
      if (length(dataFrameTotal[dataFrameTotal$id == idG & dataFrameTotal$type == typeCell & dataFrameTotal$day == dayDE ,]$log2FoldChange) > 1){
        print(dayDE)
        print(typeCell)
        print(unname(comb_size(Sup2Combi)[num_combi]))
        print(dataFrameTotal[dataFrameTotal$id == idG & dataFrameTotal$type == typeCell ,]$log2FoldChange)
        print(idG)
      }
      
    }
    num_combi=num_combi+ 1
  }
}

Violindata<-data.frame(ViolinDay=ViolinDay,
                       ViolinCell=ViolinCell,
                       ViolinName=ViolinName,
                       ViolinValue=ViolinValue)
PlotDEGonlyOneDayVioline<-ggplot(Violindata, aes(x=ViolinCell ,y=ViolinValue, fill=ViolinCell)) + # fill=name allow to automatically dedicate a color for each group
  geom_jitter(shape=16, position=position_jitter(0.4), aes(color=ViolinCell,alpha=0.001), show.legend=F)+geom_violin(scale="count",aes(alpha=0.001),show.legend = c("alpha"=F))+
  scale_fill_manual(values =c("#0072B2","#F0E442" ,"#D55E00","#CC79A7","#009E73","#56B4E9"))+scale_color_manual(values =c("#0072B2","#F0E442" ,"#D55E00","#CC79A7","#009E73","#56B4E9"))+
  coord_flip() + scale_x_discrete(limits=rev) +facet_grid(cols=vars(ViolinDay), switch = "both") +
  theme(panel.grid.major=element_blank(),axis.text.x = element_text(angle = 45)) +ylab("log2FoldChange of Gene padj <0.05")+xlab("Type Cell") +
  theme_minimal()
PlotDEGonlyOneDayVioline

save_plot(paste0(odir,analyse,"PlotDEGsigniDayTypeVioline.png"),PlotDEGonlyOneDayVioline,base_width = 10, base_height = 7)

###  CALCULATE
# Found the distribution of genes significatifs on Days Conditions by type cell
# Visualization for genes DE on only 1 days 
#
# ====================== Calculate =============================
dataFrameTotal=data.frame(signiDEgene,nbDayDE=rep(0,length(signiDEgene$baseMean)),Dynamic=rep(0,length(signiDEgene$baseMean)))
daystable<-c()
celltypetable<-c()
weightNode<-c()
ViolinDay<-c()
ViolinCell<-c()
ViolinName<-c()
ViolinValue<-c()
for ( typeCell in vectTypeCell ) {
  #keep table with signi gene
  tempoSigniDEgene<- signiDEgene %>% filter(type==typeCell)
  #keep only uniqueID sort by days
  uniqueIdByDay = lapply(unique(tempoSigniDEgene$day), function(x) filter(tempoSigniDEgene, day == x) %>% dplyr::select(id) %>% unlist() %>% as.character())
  names(uniqueIdByDay) = unique(tempoSigniDEgene$day)
  #found all combinaison where one same gene is on 1 or 2 or 3 or 4 days 
  allcombinaison =ComplexHeatmap::make_comb_mat(uniqueIdByDay)
  #keep only gene DE one only one DAY
  Sup2Combi<-allcombinaison[comb_degree(allcombinaison) == 1]
  
  #prepare vector to plot
  Node_Target<-names(comb_size(Sup2Combi))
  num_combi<-1
  for ( name in Node_Target){
    #Make table
    particulName<-str_extract_all(name,boundary("character"))
    if ( length(particulName[[1]]) <4 ){ vectDay = vectDay2} else { vectDay = vectDay1 }
    for (i in 1:(length(particulName[[1]]))){
      if ( particulName[[1]][[i]] == "1" ){
        dayDE<-vectDay[i]
        daystable<-c(daystable,vectDay[i])
        celltypetable<-c(celltypetable, typeCell)
        weightNode<-c(weightNode,unname(comb_size(Sup2Combi)[num_combi]))
      }
    }
    for (idG in extract_comb(Sup2Combi,name) ){
      dataFrameTotal[dataFrameTotal$id == idG & dataFrameTotal$type == typeCell ,]$nbDayDE<-1
      dataFrameTotal[dataFrameTotal$id == idG  & dataFrameTotal$type == typeCell ,]$Dynamic<-dayDE
      ViolinDay<-c(ViolinDay,dayDE)
      ViolinCell<-c(ViolinCell,typeCell)
      ViolinName<-c(ViolinName,unname(comb_size(Sup2Combi)[num_combi]))
      ViolinValue<-c(ViolinValue,dataFrameTotal[dataFrameTotal$id == idG & dataFrameTotal$type == typeCell ,]$log2FoldChange)
      if (length(dataFrameTotal[dataFrameTotal$id == idG & dataFrameTotal$type == typeCell ,]$log2FoldChange) > 1){
        print(dayDE)
        print(typeCell)
        print(unname(comb_size(Sup2Combi)[num_combi]))
        print(dataFrameTotal[dataFrameTotal$id == idG & dataFrameTotal$type == typeCell ,]$log2FoldChange)
      }

    }
    num_combi=num_combi+ 1
  }
}
GeneDEuniqueCond <- data.frame(
  Days=daystable,
  CellType=celltypetable,
  nbGene=weightNode
)
#cellcolors = factor( c("#0072B2","#F0E442", "#D55E00","#CC79A7","#009E73", "#56B4E9"), levels=c("ECs","FAPs","M1" ,"M2","Neutro" ,"sCs") )
PlotDEGonlyOneDay<-ggplot(GeneDEuniqueCond, aes(Days,CellType,size=nbGene,color=CellType))+geom_point()+ scale_color_manual(values =c("#0072B2","#F0E442" ,"#D55E00","#CC79A7","#009E73","#56B4E9"))  +
  scale_y_discrete(limits=rev)+   theme_ipsum()
PlotDEGonlyOneDay
save_plot(paste0(odir,analyse,"PlotDEGOnlyOneD.png"),PlotDEGonlyOneDay,base_width = 4, base_height = 7)

Violindata<-data.frame(ViolinDay=ViolinDay,
                       ViolinCell=ViolinCell,
                       ViolinName=ViolinName,
                       ViolinValue=ViolinValue)
PlotDEGonlyOneDayVioline<-ggplot(Violindata, aes(x=ViolinCell ,y=ViolinValue, fill=ViolinCell)) + # fill=name allow to automatically dedicate a color for each group
  geom_jitter(shape=16, position=position_jitter(0.4), aes(color=ViolinCell,alpha=0.001), show.legend=F)+geom_violin(scale="count",aes(alpha=0.001),show.legend = c("alpha"=F))+
  scale_fill_manual(values =c("#0072B2","#F0E442" ,"#D55E00","#CC79A7","#009E73","#56B4E9"))+scale_color_manual(values =c("#0072B2","#F0E442" ,"#D55E00","#CC79A7","#009E73","#56B4E9"))+
  coord_flip() + scale_x_discrete(limits=rev) +facet_grid(cols=vars(ViolinDay), switch = "both") +
  theme(panel.grid.major=element_blank(),axis.text.x = element_text(angle = 45)) +ylab("log2FoldChange of Gene padj <0.05")+xlab("Type Cell") +
  theme_minimal()
PlotDEGonlyOneDayVioline

save_plot(paste0(odir,analyse,"PlotDEGOnlyOneD.png"),PlotDEGonlyOneDay,base_width = 4, base_height = 7)
save_plot(paste0(odir,analyse,"PlotDEGonlyOneDayVioline.png"),PlotDEGonlyOneDay,base_width = 4, base_height = 7)


#### Same with UP & DOWN info
dataFrameTotalUPDOWN=data.frame(signiDEgene,nbDayDE=rep(0,length(signiDEgene$baseMean)),Dynamic=rep(0,length(signiDEgene$baseMean)))
daystable<-c()
celltypetable<-c()
weightNode<-c()
UPDOWN<-c()
for ( typeCell in vectTypeCell ) {
  tempoSigniDEgene<- signiDEgene %>% filter(type==typeCell)
  uniqueIdByDayUP = lapply(unique(tempoSigniDEgene$day), function(x) filter(tempoSigniDEgene, day == x)  %>% filter(log2FoldChange > 0) %>% dplyr::select(id) %>% unlist() %>% as.character())
  names(uniqueIdByDayUP) = unique(paste0(tempoSigniDEgene$day,"_UP"))
  uniqueIdByDayDOWN = lapply(unique(tempoSigniDEgene$day), function(x) filter(tempoSigniDEgene, day == x)  %>% filter(log2FoldChange < 0) %>% dplyr::select(id) %>% unlist() %>% as.character())
  names(uniqueIdByDayDOWN) = unique(paste0(tempoSigniDEgene$day,"_DOWN"))
  uniqueIdByDay<-c(uniqueIdByDayUP,uniqueIdByDayDOWN)
  #UpSetR::upset(UpSetR::fromList(uniqueIdByDay))
  allcombinaison =ComplexHeatmap::make_comb_mat(uniqueIdByDay)
  Sup2Combi<-allcombinaison[comb_degree(allcombinaison) == 1]
  Node_Target<-names(comb_size(Sup2Combi))
  num_combi<-1 #for refind number gene on this combinaison
  for ( name in Node_Target){
    particulName<-str_extract_all(name,boundary("character"))
    nbDay=length(particulName[[1]])/2
    if ( nbDay<4 ){ vectDay = vectDay2} else { vectDay = vectDay1 }
    for (i in 1:nbDay){
      if ( particulName[[1]][[i]] == "1" ){
        dayDE<-paste0(vectDay[i],"_UP")
        daystable<-c(daystable,paste0(vectDay[i],"_UP"))
        celltypetable<-c(celltypetable, typeCell)
        UPDOWN<-c(UPDOWN,"UP")
        weightNode<-c(weightNode,unname(comb_size(Sup2Combi)[num_combi]))
      }else { if (particulName[[1]][[i+nbDay]] == "1"){
        dayDE<-paste0(vectDay[i],"_DOWN")
        daystable<-c(daystable,paste0(vectDay[i],"_DOWN"))
        celltypetable<-c(celltypetable, typeCell)
        weightNode<-c(weightNode,unname(comb_size(Sup2Combi)[num_combi]))
        UPDOWN<-c(UPDOWN,"DOWN")
      }}
    } 
    num_combi=num_combi+ 1
    for (idG in extract_comb(Sup2Combi,name) ){
      dataFrameTotalUPDOWN[dataFrameTotalUPDOWN$id == idG & dataFrameTotalUPDOWN$type == typeCell ,]$nbDayDE<-1
      dataFrameTotalUPDOWN[dataFrameTotalUPDOWN$id == idG  & dataFrameTotalUPDOWN$type == typeCell ,]$Dynamic<-dayDE
    }
  }
}

GeneDEuniqueCondUPDOWN <- data.frame(
  Days=daystable,
  CellType=celltypetable,
  nbGene=weightNode,
  DiffSens=UPDOWN
)
#cellcolors = factor( c("#0072B2","#F0E442", "#D55E00","#CC79A7","#009E73", "#56B4E9"), levels=c("ECs","FAPs","M1" ,"M2","Neutro" ,"sCs") )
PlotDEGonlyOneDayUPDOWN<-ggplot(GeneDEuniqueCondUPDOWN, aes(Days,CellType,size=nbGene,color=DiffSens))+geom_point()+ scale_color_manual(values =c('red','blue'))  +
  scale_y_discrete(limits=rev)+   theme_ipsum()
PlotDEGonlyOneDayUPDOWN
save_plot(paste0(odir,analyse,"PlotDEGOnlyOneDayUPDOWN.png"),PlotDEGonlyOneDayUPDOWN,base_width = 4, base_height = 7)

###  CALCULATE
# Found the distribution of genes significatifs on Days Conditions by type cell
# Visualization with sankeyNetwork for genes DE on more 1 days
#
# ====================== Calculate =============================
sourceNode<-c()
targetNode<-c()
valueEdge<-c()
Nbconnection<-c()
for ( typeCell in vectTypeCell ) {
  tempoSigniDEgene<- signiDEgene %>% filter(type==typeCell)
  uniqueIdByDay = lapply(unique(tempoSigniDEgene$day), function(x) filter(tempoSigniDEgene, day == x) %>% dplyr::select(id) %>% unlist() %>% as.character())
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
    for (idG in extract_comb(Sup2Combi,name) ){
      dataFrameTotal[dataFrameTotal$id == idG & dataFrameTotal$type == typeCell ,]$nbDayDE<-str_count(name,"1")
      dataFrameTotal[dataFrameTotal$id == idG  & dataFrameTotal$type == typeCell ,]$Dynamic<-str_c(unique(dayDE),collapse = ".")
    }
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

my_hist<-ggplot(links,aes(value,fill = group))+geom_bar()+scale_fill_manual(values = c("white","#D3DDDC","#C7B19C","#446455"), name ="Number of Condtion DayPostInjury which Genes are DE" )
legend<-cowplot::get_legend(my_hist)
grid.arrange(legend)
saveWidget(p, file=paste0(odir, analyse,"sankeyColor3.html"))
webshot(paste0(odir,analyse, "sankeyColor3.html"), paste0(odir,analyse,"sankeyCollor3.png"), delay = 0.2)

###  CALCULATE
# Found the distribution of genes significatifs on Days Conditions by type cell
# Visualization with sankeyNetwork for genes DE on 2 conditions, 3 or 4 conditions, look afters contrast Up Down
#
# ====================== Calculate =============================
#1 Extract significative genes
#Fill dataFrameTotalUPDOWN
for ( typeCell in vectTypeCell ) {
  tempoSigniDEgene<- signiDEgene %>% filter(type==typeCell)
  uniqueIdByDayUP = lapply(unique(tempoSigniDEgene$day), function(x) filter(tempoSigniDEgene, day == x)  %>% filter(log2FoldChange > 0) %>% dplyr::select(id) %>% unlist() %>% as.character())
  names(uniqueIdByDayUP) = unique(paste0(tempoSigniDEgene$day,"_UP"))
  uniqueIdByDayDOWN = lapply(unique(tempoSigniDEgene$day), function(x) filter(tempoSigniDEgene, day == x)  %>% filter(log2FoldChange < 0) %>% dplyr::select(id) %>% unlist() %>% as.character())
  names(uniqueIdByDayDOWN) = unique(paste0(tempoSigniDEgene$day,"_DOWN"))
  uniqueIdByDay<-c(uniqueIdByDayUP,uniqueIdByDayDOWN)
  allcombinaison =ComplexHeatmap::make_comb_mat(uniqueIdByDay)
  Sup2Combi<-allcombinaison[comb_degree(allcombinaison) >= 2]
  Node_Target<-names(comb_size(Sup2Combi))
  CombNameAssocName<-set_name(Sup2Combi)
  nbDay=length(CombNameAssocName)/2
  for ( name in Node_Target){
    particulName<-str_extract_all(name,boundary("character"))
    dayDE<-c()
    if ( nbDay <4 ){ vectDay = vectDay2} else { vectDay = vectDay1 }
    if (particulName[[1]][1] == "1"){
      dayDE<-c(dayDE,paste0(vectDay[1],"_UP"))
    }else {if (particulName[[1]][1+nbDay] == "1" ){
      dayDE<-c(dayDE,paste0(vectDay[1],"_DOWN"))}}
    if (nbDay != 2){
      for(i in 2:(nbDay-1)){
        if (particulName[[1]][i] == "1"){
          dayDE<-c(dayDE,paste0(vectDay[i],"_UP"))
        }else {if (particulName[[1]][i+nbDay] == "1" ){
          dayDE<-c(dayDE,paste0(vectDay[i],"_DOWN"))}
        }}}
    else{i=1}
    if (particulName[[1]][nbDay] == "1"){
      dayDE<-c(dayDE,paste0(vectDay[nbDay],"_UP"))
    }else {if (particulName[[1]][i+nbDay+1] == "1" ){
      dayDE<-c(dayDE,paste0(vectDay[nbDay],"_DOWN"))}}
    for (idG in extract_comb(Sup2Combi,name) ){
      dataFrameTotalUPDOWN[dataFrameTotalUPDOWN$id == idG & dataFrameTotalUPDOWN$type == typeCell ,]$nbDayDE<-str_count(name,"1")
      dataFrameTotalUPDOWN[dataFrameTotalUPDOWN$id == idG  & dataFrameTotalUPDOWN$type == typeCell ,]$Dynamic<-str_c(unique(dayDE),collapse = ".")
      }
  }
}
#Prepare vector for plot + plot
dynamique <- function(tableSigni,vectTypeCell,outputprefix,oneTypeCell){
  sourceNode<-c()
  targetNode<-c()
  valueEdge<-c()
  groupetypeCell<-c()
  Nbconnection<-c()
  for ( typeCell in vectTypeCell ) {
    tempoSigniDEgene<- tableSigni %>% filter(type==typeCell)
    uniqueIdByDayUP = lapply(unique(tempoSigniDEgene$day), function(x) filter(tempoSigniDEgene, day == x)  %>% filter(log2FoldChange > 0) %>% dplyr::select(id) %>% unlist() %>% as.character())
    names(uniqueIdByDayUP) = unique(paste0(tempoSigniDEgene$day,"_UP"))
    uniqueIdByDayDOWN = lapply(unique(tempoSigniDEgene$day), function(x) filter(tempoSigniDEgene, day == x)  %>% filter(log2FoldChange < 0) %>% dplyr::select(id) %>% unlist() %>% as.character())
    names(uniqueIdByDayDOWN) = unique(paste0(tempoSigniDEgene$day,"_DOWN"))
    uniqueIdByDay<-c(uniqueIdByDayUP,uniqueIdByDayDOWN)
    allcombinaison =ComplexHeatmap::make_comb_mat(uniqueIdByDay)
    Sup2Combi<-allcombinaison[comb_degree(allcombinaison) >= 2]
    Node_Target<-names(comb_size(Sup2Combi))
    CombNameAssocName<-set_name(Sup2Combi)
    nbDay=length(CombNameAssocName)/2
    names(CombNameAssocName)<-seq(1,length(set_size(Sup2Combi)))
    num_combi<-1
    for ( name in Node_Target){
      particulName<-str_extract_all(name,boundary("character"))
      groupetypeCell<-c(groupetypeCell,rep(typeCell,nbDay-1))
      Nbconnection<-c(Nbconnection,rep(str_count(name,"1"),nbDay-1))
      if ( nbDay <4 ){ vectDay = vectDay2} else { vectDay = vectDay1 }
      if (particulName[[1]][1] == "1"){
        sourceNode<-c(sourceNode,paste0(vectDay[1],"_UP"))
      }else {if (particulName[[1]][1+nbDay] == "1" ){
        sourceNode<-c(sourceNode,paste0(vectDay[1],"_DOWN"))}
      else { 
        sourceNode<-c(sourceNode,paste0(vectDay[1],"_NODIFF"))
        }}
      if (nbDay != 2){
      for(i in 2:(nbDay-1)){
        namespl<-str_split(CombNameAssocName[i],'_')
        if (particulName[[1]][i] == "1"){
          sourceNode<-c(sourceNode,paste0(vectDay[i],"_UP"))
          targetNode<-c(targetNode,paste0(vectDay[i],"_UP"))
          valueEdge<-c(valueEdge,unname(comb_size(Sup2Combi)[num_combi]))
        }else {if (particulName[[1]][i+nbDay] == "1" ){
          sourceNode<-c(sourceNode,paste0(vectDay[i],"_DOWN"))
          targetNode<-c(targetNode,paste0(vectDay[i],"_DOWN"))
          valueEdge<-c(valueEdge,unname(comb_size(Sup2Combi)[num_combi]))}
        else { 
          sourceNode<-c(sourceNode,paste0(vectDay[i],"_NODIFF"))
          targetNode<-c(targetNode,paste0(vectDay[i],"_NODIFF"))
          valueEdge<-c(valueEdge,unname(comb_size(Sup2Combi)[num_combi]))}
        }}}
      else{i=1}
      if (particulName[[1]][nbDay] == "1"){
        targetNode<-c(targetNode,paste0(vectDay[nbDay],"_UP"))
        valueEdge<-c(valueEdge,unname(comb_size(Sup2Combi)[num_combi]))
      }else {if (particulName[[1]][i+nbDay+1] == "1" ){
        targetNode<-c(targetNode,paste0(vectDay[nbDay],"_DOWN"))
        valueEdge<-c(valueEdge,unname(comb_size(Sup2Combi)[num_combi]))}
      else { 
        targetNode<-c(targetNode,paste0(vectDay[nbDay],"_NODIFF"))
        valueEdge<-c(valueEdge,unname(comb_size(Sup2Combi)[num_combi]))}
      }
      num_combi=num_combi+ 1
    }
}
  # Make a connection data frame
  if (oneTypeCell == TRUE){
    links <- data.frame(
      source=sourceNode,
      target=targetNode,
      value=valueEdge,
      group=Nbconnection
    )
    links<-  aggregate(links$value,by=list("source"=links$source,"target"=links$target,"group"=links$group),FUN=sum)
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
    links$group <- as.factor(links$group)
    # prepare color scale: I give one specific color for each node.
    my_color <- 'd3.scaleOrdinal() .domain(["2","3","4","0",
    "D0_UP", "D2_UP","D4_UP","D7_UP", "D0_DOWN", "D2_DOWN", "D4_DOWN","D7_DOWN","D0_NODIFF", "D2_NODIFF","D4_NODIFF","D7_NODIFF",]) .range(
    ["#D3DDDC","#C7B19C","#446455","white","red","red","red","red","blue","blue","blue","blue","pink","pink","pink","pink"])'
    
  }else{
  links <- data.frame(
    source=sourceNode,
    target=targetNode,
    value=valueEdge,
    group=groupetypeCell
  )
  links<-  aggregate(links$value,by=list("source"=links$source,"target"=links$target,"group"=links$group),FUN=sum)
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
  links$group <- as.factor(links$group)
  # prepare color scale: I give one specific color for each node.
  
  my_color <- 'd3.scaleOrdinal() .domain([
  "D0_UP", "D2_UP","D4_UP","D7_UP", "D0_DOWN", "D2_DOWN", "D4_DOWN","D7_DOWN","D0_NODIFF", "D2_NODIFF","D4_NODIFF","D7_NODIFF",
  "ECs", "FAPs","sCs","M1", "M2"]) 
  .range(["red","red","red","red","blue","blue","blue","blue","pink","pink","pink","pink","#0072B2","#F0E442","#56B4E9","#D55E00","#CC79A7"])'
  # 
  # Make the Network. I call my colour scale with the colourScale argument
  
  }
  
  p <- sankeyNetwork(Links = links, Nodes = nodes, Source = "IDsource", Target = "IDtarget", 
                     Value = "x", NodeID = "name", colourScale=my_color, LinkGroup="group", NodeGroup="group", sinksRight = FALSE,
                     nodePadding=50, nodeWidth = 5,fontSize=10,iterations = 55,margin = list("right"=250), height = 650, width = 850)
  p
  saveWidget(p, file=paste0(odir, analyse,outputprefix,"sankeyColorUPDOWN.html"))
  webshot(paste0(odir,analyse,outputprefix, "sankeyColorUPDOWN.html"), paste0(odir,analyse,outputprefix,"sankeyCollorUPDOWN.png"), delay = 0.2)
  return(p)
}


palltypecell<-dynamique(signiDEgene,vectTypeCell,"AlltypeCell",FALSE)
palltypecell
psCs<-dynamique(signiDEgene,"sCs","sCs",TRUE)
psCs
pECs<-dynamique(signiDEgene,"ECs","ECs",TRUE)
pECs
pFAPs<-dynamique(signiDEgene,"FAPs","FAPs",TRUE)
pFAPs
pM2<-dynamique(signiDEgene,"M2","M2",TRUE)
pM2
pM1<-dynamique(signiDEgene,"M1","M1",TRUE)
pM1


concatMeanCountNormalized<-function(age){
  tempoTable<-read.table(paste0('data/meanCountNormalised_',age,'_D0.ECs.txt'))
  nbligne=length(tempoTable$GeneID)
  GeneID=tempoTable$GeneID
  age.meanCountNormalizedtable<-data_frame(tempoTable,type = rep("ECs",nbligne),days = rep("D0",nbligne))
  for ( t in Typecellv){
    for (day in daysv[[t]]){
      tempoTable<-read.table(paste0('data/meanCountNormalised_',age,'_',day,'.',t,'.txt'))
      nbligne=length(tempoTable$GeneID)
      if(GeneID == tempoTable$GeneID){
        tempo<-data_frame(tempoTable,type = rep(t,nbligne),days = rep(day,nbligne))
        age.meanCountNormalizedtable<-rbind(age.meanCountNormalizedtable,tempo)
      }
    }
  }
  age.meanCountNormalizedtable <- age.meanCountNormalizedtable %>% mutate( symbol = genes_df[match(age.meanCountNormalizedtable$GeneID, genes_df$Geneid),]$symbol)
  age.meanCountNormalizedtable <- age.meanCountNormalizedtable %>% unite( 'GeneDayType' , symbol:days:type ,sep = '_',remove = F)
  return(age.meanCountNormalizedtable)
}

CountNormalizedYoung<-concatMeanCountNormalized('Young')
CountNormalizedOld<-concatMeanCountNormalized('Old')

dataFrameTotalUPDOWN <- dataFrameTotalUPDOWN %>% mutate(GeneDayType=paste0(symbol,'_',day,'_',type))

dataFrameTotalUPDOWN<- dataFrameTotalUPDOWN %>% mutate( meanCountNormalizedYoung = as.numeric(CountNormalizedYoung[match(dataFrameTotalUPDOWN$GeneDayType,CountNormalizedYoung$GeneDayType),]$countNormalized)) 
dataFrameTotalUPDOWN<- dataFrameTotalUPDOWN %>% mutate( meanCountNormalizedOld = as.numeric(CountNormalizedOld[match(dataFrameTotalUPDOWN$GeneDayType,CountNormalizedOld$GeneDayType),]$countNormalized)) 
dataFrameTotalUPDOWN

saveRDS(dataFrameTotal,paste0(odir,"rds/shot_rds_full_nbDayDE.rds"))
saveRDS(dataFrameTotalUPDOWN,paste0(odir,"rds/shot_rds_full_Dynamics.rds"))


###########################################"
#Statistiques genes founds only in one days
############################################

TotUniqueGene<-Reduce(`+`,lapply(unique(fullDEsta$type), function(t) filter(fullDEsta,type == t) %>% dplyr::select(id) %>% unlist() %>% as.character() %>% unique() %>% length()))
TotUniqueGeneDE<-Reduce(`+`,lapply(unique(signiDEgene$type), function(t) filter(signiDEgene,type == t) %>% dplyr::select(id) %>% unlist() %>% as.character()%>% unique %>% length()))
#Exemple gene DE signi in more of one type cellulaire 
TotGeneDED0<-length(signiDEgene[signiDEgene$day == 'D0',]$id)
TotUniqueGeneDED0<-length(unique(signiDEgene[signiDEgene$day == 'D0',]$id))
#In which type cell there are a number max of gene DE in one days
maxDO<-GeneDEuniqueCond[GeneDEuniqueCond$nbGene == max(GeneDEuniqueCond[GeneDEuniqueCond$Days == "D0",]$nbGene),]
maxD2<-GeneDEuniqueCond[GeneDEuniqueCond$nbGene == max(GeneDEuniqueCond[GeneDEuniqueCond$Days == "D2" ,]$nbGene),]
maxD4<-GeneDEuniqueCond[GeneDEuniqueCond$nbGene == max(GeneDEuniqueCond[GeneDEuniqueCond$Days == "D4" ,]$nbGene),]
maxD7<-GeneDEuniqueCond[GeneDEuniqueCond$nbGene == max(GeneDEuniqueCond[GeneDEuniqueCond$Days == "D7",]$nbGene),]

# Number Gene DE on Only one TypeCell and the percent 
dataframeDEGonecond <- lapply(unique(dataFrameTotal$type), function(t) filter(dataFrameTotal,type == t) %>% filter(nbDayDE == 1))
names(dataframeDEGonecond)<-unique(dataFrameTotal$type)
nbUniqGeneDEonecond<-  Reduce(`+`,lapply(dataframeDEGonecond, function(d) unique(d$id) %>% length()))

dataframeDEGtwocond <- lapply(unique(dataFrameTotal$type), function(t) filter(dataFrameTotal,type == t) %>% filter(nbDayDE == 2))
names(dataframeDEGtwocond)<-unique(dataFrameTotal$type)
nbUniqGeneDEtwocond<-  Reduce(`+`,lapply(dataframeDEGtwocond, function(d) unique(d$id) %>% length()))

dataframeDEGthreecond <- lapply(unique(dataFrameTotal$type), function(t) filter(dataFrameTotal,type == t) %>% filter(nbDayDE == 3))
names(dataframeDEGthreecond)<-unique(dataFrameTotal$type)
nbUniqGeneDEthreecond<-  Reduce(`+`,lapply(dataframeDEGthreecond, function(d) unique(d$id) %>% length()))

dataframeDEGfourcond <- lapply(unique(dataFrameTotal$type), function(t) filter(dataFrameTotal,type == t) %>% filter(nbDayDE == 4))
names(dataframeDEGfourcond)<-unique(dataFrameTotal$type)
nbUniqGeneDEfourcond<-  Reduce(`+`,lapply(dataframeDEGfourcond, function(d) unique(d$id) %>% length()))

nbUniqGeneDEonecondD0 <- Reduce(`+`,lapply(dataframeDEGonecond, function(d) filter(d,day == "D0") %>% dplyr::select(id) %>% unlist() %>% as.character() %>% unique() %>% length()))
nbUniqGeneDEonecondD2 <- Reduce(`+`,lapply(dataframeDEGonecond, function(d) filter(d,day == "D2") %>% dplyr::select(id) %>% unlist() %>% as.character() %>% unique() %>% length()))
nbUniqGeneDEonecondD4 <- Reduce(`+`,lapply(dataframeDEGonecond, function(d) filter(d,day == "D4") %>% dplyr::select(id) %>% unlist() %>% as.character() %>% unique() %>% length()))
nbUniqGeneDEonecondD7 <- Reduce(`+`,lapply(dataframeDEGonecond, function(d) filter(d,day == "D7") %>% dplyr::select(id) %>% unlist() %>% as.character() %>% unique() %>% length()))
tablePercentgeneDEonecond<- data_frame("D0"=c(nbUniqGeneDEonecondD0,
                                              round(nbUniqGeneDEonecondD0/nbUniqGeneDEonecond*100, digits =0)),
                                       "D2"=c(nbUniqGeneDEonecondD2,
                                              round(nbUniqGeneDEonecondD2/nbUniqGeneDEonecond*100, digits =0)),
                                       "D4"=c(nbUniqGeneDEonecondD4,
                                              round(nbUniqGeneDEonecondD4/nbUniqGeneDEonecond*100, digits =0)),
                                       "D7"=c(nbUniqGeneDEonecondD7,
                                              round(nbUniqGeneDEonecondD7/nbUniqGeneDEonecond*100, digits =0)),
                                       "names"=c("Number Gene DE on Only one TypeCell","Percent Gene DE on Only one TypeCell")
)
tablePercentgeneDEonecond=column_to_rownames(tablePercentgeneDEonecond, var = "names") 
# Number Gene DE on 2, 3 and 4 TypeCell and the percent
tablePercentgeneDEacrosscond<- data_frame("one Day"= c(nbUniqGeneDEonecond,round(nbUniqGeneDEonecond/TotUniqueGeneDE*100,digits = 0)),
                                          "two Day"=c(nbUniqGeneDEtwocond,round(nbUniqGeneDEtwocond/TotUniqueGeneDE*100,digits = 0)),
                                          "three Day"=c(nbUniqGeneDEthreecond,round(nbUniqGeneDEthreecond/TotUniqueGeneDE*100, digits = 0)),
                                          "four Day"=c(nbUniqGeneDEfourcond, round(nbUniqGeneDEfourcond/TotUniqueGene*100, digits = 0)),
                                          "names"=c("Number Gene DE on ","Pourcent Gene DE on"))
tablePercentgeneDEacrosscond=column_to_rownames(tablePercentgeneDEacrosscond, var = "names")
#In which type cell there are a number max of gene DE in 2  or 3 or 4 days

i=0
for (t in dataframeDEGtwocond){ i=i+1; print(i) ; if (dim(t)[1] <1){ dataframeDEGtwocond<-dataframeDEGtwocond[-i] ; i=i-1 }}
i=0
for (t in dataframeDEGthreecond){ i=i+1; print(i) ; if (dim(t)[1] <1){dataframeDEGthreecond<-dataframeDEGthreecond[-i]; i=i-1}}
i=0
for (t in dataframeDEGfourcond){ i=i+1; print(i) ; if (dim(t)[1] <1){dataframeDEGfourcond<-dataframeDEGfourcond[-i]; i=i-1}}

maxGE2cond<-lapply(names(dataframeDEGtwocond), function(x) dataframeDEGtwocond[[x]] %>% dplyr::select(Dynamic) %>% unlist() %>% as.character() %>% table() %>% as_data_frame() %>% top_n(1) %>% mutate(n,"numberGene"=n/2))
names(maxGE2cond) = names(dataframeDEGtwocond)
maxGE3cond<-lapply(names(dataframeDEGthreecond), function(x) dataframeDEGthreecond[[x]] %>% dplyr::select(Dynamic) %>% unlist() %>% as.character() %>% table() %>% as_data_frame() %>% top_n(1) %>% mutate(n,"numberGene"=n/3))
names(maxGE3cond) = names(dataframeDEGthreecond)
maxGE4cond<-lapply(names(dataframeDEGfourcond), function(x) dataframeDEGfourcond[[x]] %>% dplyr::select(Dynamic) %>% unlist() %>% as.character() %>% table() %>% as_data_frame() %>% top_n(1) %>% mutate(n,"numberGene"=n/4))
names(maxGE4cond) = names(dataframeDEGfourcond)
tableMaxGE<-data_frame("ECs.Dynamics"=c(maxGE2cond$ECs$.,maxGE3cond$ECs$.,maxGE4cond$ECs$.),
                       "ECs.number"=c(maxGE2cond$ECs$numberGene,maxGE3cond$ECs$numberGene,maxGE4cond$ECs$numberGene),
                       "FAPs.Dynamics"=c(maxGE2cond$FAPs$.,maxGE3cond$FAPs$.,maxGE4cond$FAPs$.),
                       "FAPs.number"=c(maxGE2cond$FAPs$numberGene,maxGE3cond$FAPs$numberGene,maxGE4cond$FAPs$numberGene),
                       "sCs.Dynamics"=c(maxGE2cond$sCs$. ,maxGE3cond$sCs$.,maxGE4cond$sCs$.),
                       "sCs.number"=c(maxGE2cond$sCs$numberGene,maxGE3cond$sCs$numberGene,maxGE4cond$sCs$numberGene),
                       "M2.Dynamics"=c(maxGE2cond$M2$.,maxGE3cond$M2$.,NA),
                       "M2.number"=c(maxGE2cond$M2$numberGene,maxGE3cond$M2$numberGene,NA),
                       "M1.Dynamics"=c(maxGE2cond$M1$.,NA,NA),
                       "M1.number"=c(maxGE2cond$M1$numberGene,NA,NA),
)

dataframeDEGtwocondUPDOWN <-  lapply(unique(dataFrameTotalUPDOWN$type), function(t) filter(dataFrameTotalUPDOWN,type == t) %>% filter(nbDayDE == 1))
names(dataframeDEGtwocondUPDOWN)<-unique(dataFrameTotalUPDOWN$type)
dataframeDEGthreecondUPDOWN <- lapply(unique(dataFrameTotalUPDOWN$type), function(t) filter(dataFrameTotalUPDOWN,type == t) %>% filter(nbDayDE == 3))
names(dataframeDEGthreecondUPDOWN)<-unique(dataFrameTotalUPDOWN$type)
dataframeDEGfourcondUPDOWN <-  lapply(unique(dataFrameTotalUPDOWN$type), function(t) filter(dataFrameTotalUPDOWN,type == t) %>% filter(nbDayDE == 4))
names(dataframeDEGfourcondUPDOWN)<-unique(dataFrameTotalUPDOWN$type)

i=0
for (t in dataframeDEGtwocondUPDOWN){ i=i+1; print(i) ; if (dim(t)[1] <1){ dataframeDEGtwocondUPDOWN<-dataframeDEGtwocondUPDOWN[-i] ; i=i-1 }}
i=0
for (t in dataframeDEGthreecondUPDOWN){ i=i+1; print(i) ; if (dim(t)[1] <1){dataframeDEGthreecondUPDOWN<-dataframeDEGthreecondUPDOWN[-i]; i=i-1}}
i=0
for (t in dataframeDEGfourcondUPDOWN){ i=i+1; print(i) ; if (dim(t)[1] <1){dataframeDEGfourcondUPDOWN<-dataframeDEGfourcondUPDOWN[-i]; i=i-1}}


maxGE2condUPDOWN<-lapply(names(dataframeDEGtwocondUPDOWN), function(x) dataframeDEGtwocondUPDOWN[[x]] %>% dplyr::select(Dynamic) %>% unlist() %>% as.character() %>% table() %>% as_data_frame() %>% top_n(3,n) %>% arrange(desc(n)) %>% mutate(n,"numberGene"=n/2) )
names(maxGE2condUPDOWN) = names(dataframeDEGtwocondUPDOWN)
maxGE3condUPDOWN<-lapply(names(dataframeDEGthreecondUPDOWN), function(x) dataframeDEGthreecondUPDOWN[[x]] %>% dplyr::select(Dynamic) %>% unlist() %>% as.character() %>% table() %>% as_data_frame() %>% top_n(3,n) %>% arrange(desc(n))%>% mutate(n,"numberGene"=n/3) )
names(maxGE3condUPDOWN) = names(dataframeDEGthreecondUPDOWN)
maxGE4condUPDOWN<-lapply(names(dataframeDEGfourcondUPDOWN), function(x) dataframeDEGfourcondUPDOWN[[x]] %>% dplyr::select(Dynamic) %>% unlist() %>% as.character() %>% table() %>% as_data_frame() %>% top_n(3,n) %>% arrange(desc(n))%>% mutate(n,"numberGene"=n/4) )
names(maxGE4condUPDOWN) = names(dataframeDEGfourcondUPDOWN)

tableMaxGEUPDOWNECs<-data_frame( "ECs.Dynamics"=c(maxGE2condUPDOWN$ECs$.,maxGE3condUPDOWN$ECs$.,maxGE4condUPDOWN$ECs$.),
                                 "ECs.number"=c(maxGE2condUPDOWN$ECs$numberGene,maxGE3condUPDOWN$ECs$numberGene,maxGE4condUPDOWN$ECs$numberGene))
tableMaxGEUPDOWNFAPs<-data_frame("FAPs.Dynamics"=c(maxGE2condUPDOWN$FAPs$.,maxGE3condUPDOWN$FAPs$.,maxGE4condUPDOWN$FAPs$.),
                                 "FAPs.number"=c(maxGE2condUPDOWN$FAPs$numberGene,maxGE3condUPDOWN$FAPs$numberGene,maxGE4condUPDOWN$FAPs$numberGene))
tableMaxGEUPDOWNsCs<-data_frame("sCs.Dynamics"=c(maxGE2condUPDOWN$sCs$. ,maxGE3condUPDOWN$sCs$.,maxGE4condUPDOWN$sCs$.),
                                "sCs.number"=c(maxGE2condUPDOWN$sCs$numberGene,maxGE3condUPDOWN$sCs$numberGene,maxGE4condUPDOWN$sCs$numberGene))

tableMaxGEUPDOWNM2<-data_frame( "M2.Dynamics"=c(maxGE2condUPDOWN$M2$.,maxGE3condUPDOWN$M2$.),
                                "M2.number"=c(maxGE2condUPDOWN$M2$numberGene,maxGE3condUPDOWN$M2$numberGene))
tableMaxGEUPDOWNM1<-data_frame("M1.Dynamics"=c(maxGE2condUPDOWN$M1$.),
                               "M1.number"=c(maxGE2condUPDOWN$M1$numberGene))

