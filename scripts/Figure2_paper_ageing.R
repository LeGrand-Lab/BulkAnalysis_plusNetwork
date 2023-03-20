##
# Creation figure for paper
# Exploration project of 6 types cells expression during muscle regeneration in young and old cells in regeneration muscle tissue,
# johaGL 2021 + Pauline 2022
##
library(dplyr)
library(tidyverse)
library(ggplot2)
library(cowplot)
library(RColorBrewer)
library(ggthemes)
library(MASS)
library(pheatmap)
library(DESeq2)
library(ggsci) # publishing palettes
library(cowplot)
library(gridExtra)
library(ggforce)
library(extrafont)
library(colorblindr)
library(colorBlindness)
library(scales)
library(ComplexHeatmap)
library(patchwork)
library(facetious)
library(spaceNtime)
library(networkD3)
library(webshot)
library(htmlwidgets)
library(scales)
library(data.tree)
ttf_import(paths = "/home/bioinfo/R/fonts")
fonts()
loadfonts(device="postscript")
setwd("~/BulkAnalysis_plusNetwork/")

##################################
#### Path files to load or created
##################################

#Directories
ofig <- "plotsPrelim/"
inData<- "data/"
NormData <- "data/CountNormalised/"
odir <- "exam_INTER_conditions/static/"
resdir <-"PlotsDEG/"


#outputs:
metadata.rds <- "data/metadata.rds"

#Loads
fTPM<-readRDS("data/prefiltered_TPM.rds")
metadata <- readRDS(metadata.rds)
fullDEsta = readRDS(paste0(odir, "TableDEG/DEG_table_static_full.rds"))
signiDEgene <- fullDEsta %>% dplyr::filter(padj<=0.05)
dataFrameTotalUPDOWN<-readRDS(paste0(odir,"TableDEG/TableDynamicUpDownDEG.rds"))

tableGSEA<-"GSEA/TableGSEA/"
fullGSEAconcatfile<-paste0(odir,tableGSEA,"GSEA_table_static_full.csv" )
fullGSEAconcat<-read.csv(fullGSEAconcatfile,sep = " " )
GSEAsigni<-fullGSEAconcat %>% filter(padj<=0.05) %>% mutate(sens=ifelse(NES >0,"Up","Down"))


genes_df <- read.table("data/genesinfo.csv", sep="\t", header=T)
daysv<-list()
daysv<-lapply( unique(fullDEsta$type), function(t) unique(fullDEsta[fullDEsta$type == t,]$day) )
names(daysv)<-unique(fullDEsta$type)
Typecellv<-unique(fullDEsta$type)
days<-unique(fullDEsta$day)

##
# Correlation matrice , spearman test
##
# =========================================================================
orderTypecell=c("ECs","FAPs","MuSCs","Neutrophils","Inflammatory-Mac","Resolving-Mac")
colorsType=c("#10b387ff","#3d85c6ff","#b171f1ff","#f0e442ff","#ff9900ff","#cc0000ff")
names(colorsType) = orderTypecell
myannot = data.frame(CellType=factor(metadata$type, levels = orderTypecell ),
                           Time = factor(metadata$time,levels= str_sort(unique(metadata$time))),
                     Age = factor(metadata$age, levels = str_sort(unique(metadata$age), decreasing =T)))
rownames(myannot) = rownames(metadata)



colorsTime = c(brewer.pal(9,"BuPu")[3],brewer.pal(9,"BuPu")[5],brewer.pal(9,"BuPu")[7],brewer.pal(9,"BuPu")[9])
names(colorsTime) = sort(unique(myannot$Time))
colorsAge = c("#ffc35d","#019190")
names(colorsAge) = c("Young","Old")
mycolors = list("CellType"=colorsType,"Time" = colorsTime, "Age" =colorsAge )


lgdType<-Legend(at=levels(myannot$CellType),legend_gp = gpar(fill = colorsType), title="CellType",labels_gp = gpar (fontsize= 5),title_gp = gpar(fontsize = 6, fontface = "bold"))
lgdAge<-Legend(at=levels(myannot$Age),legend_gp = gpar(fill = colorsAge), title="Age",labels_gp = gpar (fontsize= 5),title_gp = gpar(fontsize = 6, fontface = "bold"))
lgdTime<-Legend(at=levels(myannot$Time),legend_gp = gpar(fill = colorsTime ), title="Time",labels_gp = gpar (fontsize= 5),title_gp = gpar(fontsize = 6, fontface = "bold"))

tableGSEA<-"GSEA/TableGSEA/"
fullGSEAconcatfile<-paste0(odir,tableGSEA,"GSEA_table_static_full.csv" )
fullGSEAconcat<-read.csv(fullGSEAconcatfile,sep = " " )
GSEAsigni<-fullGSEAconcat %>% filter(padj<=0.05) %>% mutate(sens=ifelse(NES >0,"Up","Down"))

#png(paste0(ofig,"Figure1_corpus.png"),units = "in", width=8, height= 11, res = 300, family = "Arial")
#ggpubr::ggarrange(CovHeatmap_grob,ggpubr::ggarrange(facet1_ld , facet2.5_ld , facet4_ld, ncol=1,heights = c(0.25,0.25,0.5), labels=c("B","C","D")),g5, facet6_ld,ncol=2,nrow = 2, labels = c("A","","E","F"))
#dev.off()

#Found genes DE on one, or more days
#Observe the dynamic of expression for Old in comparison to Young across days on one Type cell -> sankeyNetwork
#Violine plot of DEG 
signiDEgene<-signiDEgene %>% mutate(day.type=paste0(day,'_',type))
nbGeneSignificant<-table(signiDEgene$day.type)
Violindata<-data.frame(ViolinDT=paste0(signiDEgene$day,"_",signiDEgene$type),
                       ViolinDay=factor(signiDEgene$day,levels= levels(myannot$Time) ),
                       CellType=factor(signiDEgene$type, levels= levels(myannot$CellType)),
                       ViolinValue=signiDEgene$log2FoldChange,
                       nbGeneSignificant = lapply(signiDEgene$day.type, function(dt) nbGeneSignificant[[dt]]) %>% unlist(),
                       a=rep("a",length(signiDEgene$day.type)))
PlotDEGonlyOneDayVioline<-ggplot(Violindata, aes(x=a ,y=ViolinValue, fill=nbGeneSignificant)) +# fill=name allow to automatically dedicate a color for each group
  geom_jitter( shape=16, size=1, position=position_jitter(0.4), aes(color=nbGeneSignificant,alpha=0.001), show.legend=F)+
  geom_violin(scale="count",aes(alpha=0.001),show.legend = c("alpha"=F))+
  scale_fill_gradient(low = "#E7E1EF",high="#67001F",labels = scales::label_comma())+scale_color_gradient(low = "#E7E1EF",high="#67001F",labels = scales::label_comma())+
  coord_flip()+
  scale_x_discrete(limits=rev)  +   facet_grid_blank(vars(CellType),vars(ViolinDay), drop = FALSE) +
  geom_hline(yintercept = 0, color="black", linetype="dashed",lwd=0.5)+
  ylab("log2FoldChange")+xlab("Type Cell") +
  theme(axis.title = element_text(size=5, face = "bold"),axis.text.x=element_text(size=3,colour = "black", face="bold"),axis.text.y=element_text(colour = "white"),legend.title=element_text( size=6, face = "bold"),legend.text=element_text(size=4),title=element_text(size=7, face = "bold"),
        legend.spacing= unit(0.1,"points"),legend.position = "bottom",legend.spacing.y=unit(0.01,"mm"),legend.key.size = unit(3, "mm"),panel.grid =element_line(color="white"),strip.text.x = element_text(size=10,color = "white", face= "bold"),strip.text.y = element_text(size=5.5,color = "black", face= "bold"),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.background = element_blank())+labs(fill = "Count DEG padj<0.05")+
  ggtitle("Distribution of significatif differentially expressed genes Old vs Young \n according to log2Foldchange")



PlotDEGonlyOneDayVioline
g5 <- ggplot_gtable(ggplot_build(PlotDEGonlyOneDayVioline))
stripRowName <- which(grepl('strip', g5$layout$name))

k <- 1
fills <- c(colorsTime,colorsType)
for (i in stripRowName) {
  j <- which(grepl('rect', g5$grobs[[i]]$grobs[[1]]$childrenOrder))
  g5$grobs[[i]]$grobs[[1]]$children[[j]]$gp$fill <- fills[k]
  k <- k+1
}
grid::grid.draw(g5)
facet5<-grid.arrange(g5)
save_plot(paste0(ofig,"PlotDEGsigniDayTypeVioline.png"),g5,base_width = 10, base_height = 7)

### Place sankey plot
colorCombiDayDEG<-c("#fef6be","#e2f4df","#dbeaf5","#d8d9ea",
                    "#d6efa6","#c2d9ed","#ffcce4","#a0dbb5","#d5b7d7","#8895c2",
                    "#008349","#065fa3","#d42857","#7e2175","#67081d")
names(colorCombiDayDEG)<-c("D0","D2","D4","D7"
                           ,"D0.D2","D0.D4","D0.D7","D2.D4","D2.D7","D4.D7",
                           "D0.D2.D4","D0.D2.D7","D0.D4.D7","D2.D4.D7","D0.D2.D4.D7")
colorCombiDayDEG<-c("#ffe119","#bfef45","#42d4f4","#f032e6","#3cb44b","#4363d8","#911eb4","#006600","#000075","#67081d")
names(colorCombiDayDEG)<-c("D0","D2","D4","D7"
                           ,"D0.D2","D2.D4","D4.D7",
                           "D0.D2.D4","D2.D4.D7","D0.D2.D4.D7")

linksToRemove<-c("D0.D4","D0.D7","D2.D7","D0.D2.D7","D0.D4.D7")
colorsConditions<-colorCombiDayDEG[! names(colorCombiDayDEG) %in% linksToRemove]
lgddaysDEG<-Legend(at=names(colorCombiDayDEG),legend_gp = gpar(fill = colorCombiDayDEG), title="days_DEG",labels_gp = gpar (fontsize= 5),title_gp = gpar(fontsize = 6, fontface = "bold"))
lgddaysDEG_grob=grid.grabExpr(draw(lgddaysDEG)) 
grid.arrange(lgddaysDEG_grob)
grid.draw(lgddaysDEG_grob)

Violindata2<-data.frame(ViolinDay=factor(signiDEgene$day,levels= levels(myannot$Time) ),
                        CellType=factor(signiDEgene$type, levels= levels(myannot$CellType)),
                        ViolinValue=rep("b",length(signiDEgene$day.type)),
                        a=rep("a",length(signiDEgene$day.type)))
PlotDEGonlyOneDayVioline2<-ggplot(Violindata2, aes(x=a ,y=ViolinValue)) +# fill=name allow to automatically dedicate a color for each group
  geom_jitter(color="white",fill="white",show.legend=F)+
  coord_flip()+
  facet_grid_blank(vars(CellType),vars(ViolinDay), drop = FALSE) +
  ylab("Day")+xlab("Type Cell") +
  ggtitle("Flow of significatif differentially expressed gene \n on consecutive days and on one cell type")+
  theme(axis.title = element_text(size=5, face = "bold"),axis.text=element_text(colour = "white"),title=element_text(size=7, face = "bold"),
        legend.spacing= unit(0.1,"points"),legend.position = "bottom",legend.spacing.y=unit(0.01,"mm"),legend.key.size = unit(3, "mm"),panel.grid =element_line(color="white"),strip.text.x = element_text(size=10,color = "white", face= "bold"),strip.text.y = element_text(size=5.5,color = "black", face= "bold"),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.background = element_blank())



PlotDEGonlyOneDayVioline2
g6 <- ggplot_gtable(ggplot_build(PlotDEGonlyOneDayVioline2))
stripRowName <- which(grepl('strip', g6$layout$name))

k <- 1
fills <- c(colorsTime,colorsType)
for (i in stripRowName) {
  j <- which(grepl('rect', g6$grobs[[i]]$grobs[[1]]$childrenOrder))
  g6$grobs[[i]]$grobs[[1]]$children[[j]]$gp$fill <- fills[k]
  k <- k+1
}
grid::grid.draw(g6)
facet6_ld<-grid.arrange(blank,lgddaysDEG_grob,ncol=2,widths=c(7,1.5))

png(paste0(ofig,"Figure1_corpus.png"),units = "in", width=8, height= 11, res = 300, family = "Arial")
ggpubr::ggarrange(CovHeatmap_grob,ggpubr::ggarrange(facet1_ld , facet2.5_ld , facet4_ld, ncol=1,heights = c(0.25,0.25,0.5), labels=c("B","C","D")),g5, facet6_ld,ncol=2,nrow = 2, labels = c("A","","E","F"))
dev.off()


### SankeyPlot

sourceNode<-c()
targetNode<-c()
valueEdge<-c()
Nbconnection<-c()
CombiDayDEG<-c()
for ( typeCell in Typecellv ) {
  tempoSigniDEgene<- signiDEgene %>% filter(type==typeCell)
  #Get list DEG by day
  uniqueIdByDay = lapply(unique(tempoSigniDEgene$day), function(x) filter(tempoSigniDEgene, day == x) %>% dplyr::select(id) %>% unlist() %>% as.character())
  names(uniqueIdByDay) = unique(tempoSigniDEgene$day)
  
  #Get all combination DEG 
  allcombinaison =ComplexHeatmap::make_comb_mat(uniqueIdByDay)
  #Extract pattern combination
  Sup2Combi<-allcombinaison[comb_degree(allcombinaison) >= 1]
  Node_Target<-names(comb_size(Sup2Combi))
  
  num_combi<-1
  
  #For one combination
  #Extract the number of day DEG
  #Found the name of combonation with day
  for ( name in Node_Target){
    particulName<-str_extract_all(name,boundary("character"))[[1]]
    names(particulName)<-1:length(particulName)
    
    Nodes<-lapply(1:length(particulName), function(i) if(particulName[i] == "1"){daysv[[typeCell]][i]} ) %>% unlist()
    if (str_count(name,"1")==1){CombiDayDEG<-c(CombiDayDEG,rep(str_c(Nodes, collapse = "."),4))} else {CombiDayDEG<-c(CombiDayDEG,rep(CombiDayDEG<-c(str_c(Nodes, collapse = ".")),str_count(name,"1")-1))}
    if (str_count(name,"1")==1){Nbconnection<-c(Nbconnection,rep(1,4))} else {Nbconnection<-c(Nbconnection,rep(str_count(name,"1"),str_count(name,"1")-1))}
    if (length(Nodes) == 1){
      sourceNode<-c(sourceNode,paste0(Nodes[1],"_",typeCell))
      targetNode<-c(targetNode,paste0(Nodes[1],"_",typeCell,"_only1"))
      valueEdge<-c(valueEdge,unname(comb_size(Sup2Combi)[num_combi]))
      sourceNode<-c(sourceNode,paste0(Nodes[1],"_",typeCell,"_only1"))
      targetNode<-c(targetNode,paste0(Nodes[1],"_",typeCell,"_only2"))
      valueEdge<-c(valueEdge,unname(comb_size(Sup2Combi)[num_combi]))
      sourceNode<-c(sourceNode,paste0(Nodes[1],"_",typeCell,"_only2"))
      targetNode<-c(targetNode,paste0(Nodes[1],"_",typeCell,"_only3"))
      valueEdge<-c(valueEdge,unname(comb_size(Sup2Combi)[num_combi]))
      sourceNode<-c(sourceNode,paste0(Nodes[1],"_",typeCell,"_only3"))
      targetNode<-c(targetNode,paste0(Nodes[1],"_",typeCell))
      valueEdge<-c(valueEdge,unname(comb_size(Sup2Combi)[num_combi]))
    }else{
      sourceNode<-c(sourceNode,paste0(Nodes[1],"_",typeCell))
      if(length(Nodes)-1 > 1){
        for (i in 2:(length(Nodes)-1)){
          sourceNode<-c(sourceNode,paste0(Nodes[i],"_",typeCell))
          targetNode<-c(targetNode,paste0(Nodes[i],"_",typeCell))
          valueEdge<-c(valueEdge,unname(comb_size(Sup2Combi)[num_combi]))
        }}
      targetNode<-c(targetNode,paste0(Nodes[length(Nodes)],"_",typeCell))
      valueEdge <-c(valueEdge,unname(comb_size(Sup2Combi)[num_combi]))
    }
    num_combi=num_combi+ 1
  } 
}  
#color
colorNbconnection<-c("white","#FABDA6","#D3DDDC","#C7B19C","#446455")
names(colorNbconnection)<-c("0","1","2","3","4")


# Make a connection data frame
# Total links
links <- data.frame(
  from=str_replace_all(sourceNode,"-",""),
  to=str_replace_all(targetNode,"-",""),
  Nbconnection=Nbconnection,
  CombiDayDEG=CombiDayDEG,
  substance= CombiDayDEG,
  quantity=sqrt(valueEdge)
)
links<-links %>% arrange(substance)
# From these flows we need to create a node data frame: it lists every entities involved in the flow
nodes <- data.frame(
  ID=c(as.character(links$from), as.character(links$to)) %>% 
    unique() %>% str_replace_all("-","")
)
nodes$label <- as.factor(c(as.character(links$from), as.character(links$to)) %>% 
                           unique())
nodes<-nodes  %>% mutate(day=sapply(strsplit(nodes$ID,"_"), `[`, 1)) %>% mutate(type=factor(sapply(strsplit(nodes$ID,"_"), `[`, 2),levels=str_replace_all(orderTypecell,"-",""),)) %>% mutate(loops=sapply(strsplit(nodes$ID,"_"), `[`, 3))
nbDay<-length(unique(nodes$day))
nbtype<-length(unique(nodes$type))
nodes<-nodes  %>% arrange(loops,type,day)
nodes<- nodes %>% mutate(label=lapply(1:length(nodes$day) ,function(i) ifelse( is.na(nodes[i,]$loops) == T,nodes[i,]$ID,"" )) %>% unlist()) 
nodes<- nodes %>%  mutate(ID=lapply( 1:length(nodes$day) ,function(i) ifelse( is.na(nodes[i,]$loops) == F,paste0(".",nodes[i,]$ID), nodes[i,]$ID)) %>% unlist())
poslabel=c("left","above","below","right")
nodes$label_pos<-c(rep("none",54),rep("below",18))
xpos4<-c(-7,-2,3,8)
xpos3<-c(-2,3,8)
xpos2<-c(-2,3)
xpos1<-c(-2)
nodes$x<-c(rep(xpos4+2,3),xpos1+2,xpos2+2,xpos3+2,rep(xpos4,3),xpos1,xpos2,xpos3,rep(xpos4-2,3),xpos1-2,xpos2-2,xpos3-2,rep(xpos4,3),xpos1,xpos2,xpos3)
nodes<-nodes  %>% arrange(loops,day,type)
yposD2<-c(9,4,-1,-6,-11,-16)
yposD4<-c(9,4,-1,-11,-16)
yposD7<-c(9,4,-1,-16)
yposD0<-c(9,4,-1)
nodes$y<-c(yposD0+1,yposD2+1,yposD4+1,yposD7+1,yposD0+2,yposD2+2,yposD4+2,yposD7+2,yposD0+1,yposD2+1,yposD4+1,yposD7+1,yposD0,yposD2,yposD4,yposD7)
nodes$dir<-c(rep("up",18),rep("left",18),rep("down",18),rep("right",18))

dblue<-"#67001F"
# node style
ns <- list(type="arrow",gp=gpar(fill=dblue, col="white", lwd=3),
           length=0.5,
           label_gp=gpar(col="black", fontsize=6,fontface="bold"),
           mag_pos="label", mag_fmt="%.0f", mag_gp=gpar(fontsize=1,col="white"))


palettesCombiDayDEG<-data.frame(substance=names(colorCombiDayDEG),
                                color=colorCombiDayDEG)
my_title_CombiDayDEG <- paste0("")
attr(my_title_CombiDayDEG, "gp") <- grid::gpar(fontsize=12, fontface="bold", col="black")
links$substance<-links$CombiDayDEG
rowtoRemove<-subset(links, substance %in% linksToRemove) %>% rownames()%>% as.numeric()
links2=links[-rowtoRemove,]

testsankeysave<-PantaRhei::sankey(nodes, links2, palettesCombiDayDEG,
                  max_width=0.2, rmin=0.5,
                  node_style=ns,
                  page_margin=c(0, 0.05, 0, 0.05),
                  title=my_title_CombiDayDEG, legend=NULL)

png(paste0(ofig,"Figure2_tail.png"),units = "in", width=4, height= 5.5, res = 300, family = "Arial")
PantaRhei::sankey(nodes, links2, palettesCombiDayDEG,
                  max_width=0.15, rmin=0.5,
                  node_style=ns,
                  page_margin=c(0, 0.05, 0, 0.05),
                  title=my_title_CombiDayDEG, legend=NULL)
dev.off()


tableSankeyplot<- function(signiDEgene,vectTypeCell,minCombi,NODIFF){
  sourceNode<-c()
  targetNode<-c()
  valueEdge<-c()
  NbUpDown<-c()
  Nbconnection<-c()
  CombiDayDEG<-c()
  
  for ( typeCell in vectTypeCell ) {
    print(typeCell)
    #prepare list gene_SENS DEG by day
    tempoSigniDEgene<- signiDEgene %>% filter(type==typeCell)
    uniqueIdByDayUP = lapply(unique(tempoSigniDEgene$day), function(x) filter(tempoSigniDEgene, day == x)  %>% filter(log2FoldChange > 0) %>% dplyr::select(id) %>% unlist() %>% as.character())
    names(uniqueIdByDayUP) = unique(paste0(tempoSigniDEgene$day,"_UP"))
    uniqueIdByDayDOWN = lapply(unique(tempoSigniDEgene$day), function(x) filter(tempoSigniDEgene, day == x)  %>% filter(log2FoldChange < 0) %>% dplyr::select(id) %>% unlist() %>% as.character())
    names(uniqueIdByDayDOWN) = unique(paste0(tempoSigniDEgene$day,"_DOWN"))
    uniqueIdByDay<-c(uniqueIdByDayUP,uniqueIdByDayDOWN)
    
    #Get combinaison
    allcombinaison =ComplexHeatmap::make_comb_mat(uniqueIdByDay)
    Sup2Combi<-allcombinaison[comb_degree(allcombinaison) >= minCombi]
    
    #Get name combinaison
    Node_Target<-names(comb_size(Sup2Combi))
    if(length(Node_Target) < minCombi){next}
    num_combi<-1
    
    #For one combination
    #Extract the number of day DEG
    #Found the name of combination with day+sens(UP,DOWN)
    #Extract Day + sens (UP, DOWN,UNDIFF)
    for ( name in Node_Target){
      print(name)
      particulName<-str_extract_all(name,boundary("character"))[[1]]
      names(particulName)<-1:length(particulName)
      #extract day DEG UP
      tempoName1<-lapply(1:(length(particulName)/2), function(i) if(particulName[[i]] == "1"){paste0(daysv[[typeCell]][i],"_UP")} ) %>% unlist()
      #extract day DEG DOWN
      tempoName2<-lapply(((length(particulName)/2)+1):length(particulName), function(i) if(particulName[[i]] == "1"){paste0(daysv[[typeCell]][i-(length(particulName)/2)],"_DOWN")} ) %>% unlist()
      #Extract day UNDIFF
      tempoName3<-lapply(1:(length(particulName)/2), function(i) if(particulName[[i]] == "0" && particulName[[i+(length(particulName)/2)]] == "0"){paste0(daysv[[typeCell]][i],"_UNDIFF")} ) %>% unlist()
      
      #Extract only day DEG
      tempoName1bis<-lapply(1:(length(particulName)/2), function(i) if(particulName[[i]] == "1"){daysv[[typeCell]][i]} ) %>% unlist()
      tempoName2bis<-lapply(((length(particulName)/2)+1):length(particulName), function(i) if(particulName[[i]] == "1"){daysv[[typeCell]][i-(length(particulName)/2)]} ) %>% unlist()
      tempoName<-str_c(sort(c(tempoName1bis,tempoName2bis)),collapse = ".")
      
      if (NODIFF == T ){
        Nodes<-sort(c(tempoName1,tempoName2,tempoName3))
        CombiDayDEG<-c(CombiDayDEG,rep(tempoName,length(Nodes)-1))
        Nbconnection<-c(Nbconnection,rep(length(c(tempoName1,tempoName2)),length(Nodes)-1))
        NbUpDown<-c(NbUpDown,rep(paste0(ifelse(length(tempoName1)==0,"",paste0(length(tempoName1),"_UP.")),ifelse(length(tempoName2)==0,"",paste0(length(tempoName2),"_DOWN"))),length(Nodes)-1))
        
      }else{    
        Nodes<-sort(c(tempoName1,tempoName2))
        if (length(Nodes)==1){Nbconnection<-c(Nbconnection,rep(length(c(tempoName1,tempoName2)),4))} else {Nbconnection<-c(Nbconnection,rep(length(c(tempoName1,tempoName2)),length(Nodes)-1))}
        if (length(Nodes)==1){CombiDayDEG<-c(CombiDayDEG,rep(tempoName,4))} else {CombiDayDEG<-c(CombiDayDEG,rep(tempoName,length(Nodes)-1))}
        if (length(Nodes)==1){NbUpDown<-c(NbUpDown,rep(paste0(ifelse(length(tempoName1)==0,"",paste0(length(tempoName1),"_UP.")),ifelse(length(tempoName2)==0,"",paste0(length(tempoName2),"_DOWN"))),4))} else {NbUpDown<-c(NbUpDown,rep(paste0(ifelse(length(tempoName1)==0,"",paste0(length(tempoName1),"_UP.")),ifelse(length(tempoName2)==0,"",paste0(length(tempoName2),"_DOWN"))),length(Nodes)-1))}
        
      }
      
      #Fill source node, target node (day consecutive) + nb gene with this combination
      if (length(Nodes) == 1){
        sourceNode<-c(sourceNode,Nodes[1])
        targetNode<-c(targetNode,paste0(Nodes[1],"_only1"))
        valueEdge<-c(valueEdge,unname(comb_size(Sup2Combi)[num_combi]))
        sourceNode<-c(sourceNode,paste0(Nodes[1],"_only1"))
        targetNode<-c(targetNode,paste0(Nodes[1],"_only2"))
        valueEdge<-c(valueEdge,unname(comb_size(Sup2Combi)[num_combi]))
        sourceNode<-c(sourceNode,paste0(Nodes[1],"_only2"))
        targetNode<-c(targetNode,paste0(Nodes[1],"_only3"))
        valueEdge<-c(valueEdge,unname(comb_size(Sup2Combi)[num_combi]))
        sourceNode<-c(sourceNode,paste0(Nodes[1],"_only3"))
        targetNode<-c(targetNode,Nodes[1])
        valueEdge<-c(valueEdge,unname(comb_size(Sup2Combi)[num_combi]))
      }else{
        sourceNode<-c(sourceNode,Nodes[1])
        if(length(Nodes)-1 > 1){
          for (i in 2:(length(Nodes)-1)){
            sourceNode<-c(sourceNode,Nodes[i])
            targetNode<-c(targetNode,Nodes[i])
            valueEdge<-c(valueEdge,unname(comb_size(Sup2Combi)[num_combi]))
          }}
        targetNode<-c(targetNode,Nodes[length(Nodes)])
        valueEdge <-c(valueEdge,unname(comb_size(Sup2Combi)[num_combi]))
      }
      
      num_combi=num_combi+ 1
    }
    
  }
  
  links <- data.frame(
    from=str_replace_all(sourceNode,"-",""),
    to=str_replace_all(targetNode,"-",""),
    NbUpDown= NbUpDown,
    Nbconnection=Nbconnection,
    CombiDayDEG=CombiDayDEG,
    quantity=valueEdge
  )
  
  
  return(links)
}

# sankeyPantaRhei , made a sankey plot, where you can choice coordinate of nodes, label nodes,group color links by group (NbUpDown,Nbconnection,CombiDayDEG)
#links1 = output of tableSankeyplot with at least nodes from, nodes to, value of the flow, and group
#dblue= color of nodes
#group=
#NbUpDown: group links by pattern number of UP and number od DOWN
#Nbconnection: group by number of day DEG
#CombiDayDEG: group by pattern of days DEG
#colorvect= vector color with names unique(links1$group)
#outputprefix= start of plot title
#sqrt= T/F, value in square root
#NODIFF= T/F , nodes UNDIFF?
sankeyPantaRhei<-function(links1,dblue,group,colorvect,outputprefix,sqrt,NODIFF){
  if ( sqrt == T){
    links1$quantity<-sqrt(links1$quantity)
    titlesqrt<-"(sqrt)"}else{titlesqrt<-""
    }
  print(links1)
  links <- data.frame(
    from=links1$from,
    to=links1$to,
    substance= as.factor(links1 %>% dplyr::select(group) %>% unlist()),
    #substance= as.factor(links1$CombiDayDEG),
    quantity=links1$quantity
  )
  print(links)
  links<-links %>% arrange(substance)
  # From these flows we need to create a node data frame: it lists every entities involved in the flow
  nodes <- data.frame(
    ID=c(as.character(links$from), as.character(links$to)) %>% 
      unique() %>% str_replace_all("-","")
  )
  nodes$label <- as.factor(c(as.character(links$from), as.character(links$to)) %>% 
                             unique())
  nodes<-nodes  %>% mutate(day=sapply(strsplit(nodes$ID,"_"), `[`, 1)) %>% mutate(sens=sapply(strsplit(nodes$ID,"_"), `[`, 2)) %>% mutate(loops=sapply(strsplit(nodes$ID,"_"), `[`, 3))
  nbDay<-length(unique(nodes$day))
  print(nbDay)
  nodes<-nodes  %>% arrange(loops,sens,day)
  print(nodes)
  if("only1" %in% nodes$loops){only=T}else{only=F}
  if(length(nodes %>% filter(sens== "UNDIFF") %>% dplyr::select(sens) %>% unlist()) >=1){
    replabel=3
    titleNODIFF="UNDIFF"}else{
      NODIFF = F
      replabel=2
      titleNODIFF=""}
  if (nbDay == 4){
    poslabel=c("left","above","below","right")
    xpos<-c(-7,-2,3,8)
    ypos<-c(-2,-3,-4,-5)
  } else if(nbDay == 3){
    poslabel=c("left","above","right")
    xpos<-c(-7,-2,3)
    ypos<-c(-2,-3,-4)
  } else{
    poslabel=c("left","right")
    xpos<-c(-7,3)
    ypos<-c(-2,-3)
  }
  
  nodes<- nodes %>% mutate(label=lapply(1:length(nodes$sens) ,function(i) ifelse( is.na(nodes[i,]$loops) == T,nodes[i,]$ID,"" )) %>% unlist()) 
  nodes<- nodes %>%  mutate(ID=lapply( 1:length(nodes$sens) ,function(i) ifelse( is.na(nodes[i,]$loops) == F,paste0(".",nodes[i,]$ID), nodes[i,]$ID)) %>% unlist())
  if(only == T){nodes$label_pos<-c(rep("none",nbDay*6),rep(poslabel,replabel))}else{nodes$label_pos<-c(rep(poslabel,replabel))}
  nodes$label_align<-rep("",length(nodes$ID))
  if(NODIFF == T){nodes$x<-c(xpos,xpos+1,xpos)}else{ if(only == T){nodes$x<-c(xpos+2,xpos+2,xpos,xpos,xpos-2,xpos-2,rep(xpos,2))} else{nodes$x<-c(xpos,xpos)}}
  if(NODIFF == T){nodes$y<-c(ypos,rep(0,nbDay),ypos*-1)}else{ if(only == T){nodes$y<-c(ypos-1,(ypos*-1)+1,ypos-2,(ypos*-1)+2,ypos-1,(ypos*-1)+1,ypos,ypos*-1)}else{nodes$y<-c(ypos,ypos*-1)}}
  if(only == T){nodes$dir<-c(rep("down",nbDay),rep("up",nbDay),rep("left",nbDay*2),rep("up",nbDay),rep("down",nbDay),rep("right",nbDay*2))}else{nodes$dir<-rep("right",nbDay*replabel)}
  
  palettes<-data.frame(substance=names(colorvect),
                       color=colorvect)
  my_title <- paste0("")
  attr(my_title, "gp") <- grid::gpar(fontsize=6, fontface="bold", col="black")
  
  # node style
  ns <- list(type="arrow",gp=gpar(fill=dblue, col="white", lwd=3),
             length=0.5,
             label_gp=gpar(col=dblue, fontsize=6,fontface="bold"),
             mag_pos="label", mag_fmt="%.0f", mag_gp=gpar(fontsize=6,fontface="bold",col=dblue))
  
  #pdf(paste0(analyse,"test_group_",group,"_",titleNODIFF,'_',sqrt,'_',outputprefix,"_sankeyColorUPDOWN.pdf"), width=10, height=7) # Set up PDF device
  print(nodes)
  PantaRhei::sankey(nodes, links, palettes,
                    max_width=0.15, rmin=0.5,
                    node_style=ns,
                    page_margin=c(0, 0.05, 0, 0.05),
                    title=my_title, legend=F )
  #dev.off()
  
}

colorNbUpDown<-c("plum","plum4",
                 "mediumpurple1","mediumpurple4",
                 "violetred1","violetred4",
                 "LightCyan","Aquamarine","skyblue2","royalblue3",
                 "LightPink","indianred1","firebrick1","red3" )

names(colorNbUpDown)<-c("1_UP.1_DOWN" ,"2_UP.2_DOWN",
                        "1_UP.2_DOWN", "1_UP.3_DOWN" ,
                        "2_UP.1_DOWN","3_UP.1_DOWN",
                        "1_DOWN","2_DOWN","3_DOWN" ,"4_DOWN",
                        "1_UP.","2_UP.","3_UP.","4_UP.")
linksToRemove<-c("1_UP.1_DOWN","2_UP.2_DOWN","1_UP.3_DOWN","2_UP.1_DOWN","3_UP.1_DOWN","4_DOWN","4_UP.")
colorNbUpDown<-colorNbUpDown[! names(colorNbUpDown) %in% linksToRemove]
lgddaysDEGUpDown<-Legend(at=names(colorNbUpDown),legend_gp = gpar(fill = colorNbUpDown), title="Number Day which DEG \n are UP or DOWN",labels_gp = gpar (fontsize= 5),title_gp = gpar(fontsize = 6, fontface = "bold"), direction = "horizontal",nrow = 2,title_position ="leftcenter")
lgddaysDEG_UpDown_grob=grid.grabExpr(draw(lgddaysDEGUpDown)) 
grid.arrange(lgddaysDEG_UpDown_grob)
blank<-grob()
facet5.5_ld<-grid.arrange(blank,lgddaysDEG_UpDown_grob, nrow=2,heights=c(7,1.5))

tableSankeyplotFAPs<-tableSankeyplot(signiDEgene,"Resolving-Mac",1,F)

png(paste0(ofig,"SankeyResolvingMac.png"),units = "in", width=3, height=5, res = 800, family = "Arial")
sankeyPantaRhei(tableSankeyplotFAPs,colorsType[["Resolving-Mac"]],"NbUpDown",colorNbUpDown,paste0("Resolving-Mac","_combisupp1"),T,F)
dev.off()

signiDEgeneSelenop<-dataFrameTotalUPDOWN %>% filter(symbol=="Selenop",type=="Resolving-Mac") %>% pivot_longer(cols = starts_with("meanCountNormalized"),names_to = "Age",names_prefix = "meanCountNormalized",values_to ="meanCountNormalized") %>% as.data.frame() %>% mutate(Age=factor(Age,levels=c('Young','Old')))

PlotGenelogFacrossDay2 <- ggplot(signiDEgeneSelenop, aes(x=day, y=meanCountNormalized,group=paste(symbol,Age),shape=Age) )+
  geom_line() +
  scale_color_gradient2(midpoint=0,  low="blue", mid="darkgrey",high="red",name="log2FoldChange")+
  geom_hline(yintercept = 0, color="black", linetype="dashed",lwd=0.5,slope=F)+
  geom_point(aes(col=log2FoldChange,size=meanCountNormalized),show.legend = c(col=T,size=F))+
  ylab("Mean count normalised ")+ggtitle("Selenoprotein expression in resolving macrophages : Old vs Young")+
  theme_bw()+ theme(axis.title = element_text(size=5, face = "bold"),axis.text=element_text(size=4),legend.title=element_text(size=6, face = "bold"),legend.text=element_text(size=4),title=element_text(size=7, face = "bold"),
                    legend.spacing= unit(0.1,"points"),legend.spacing.y=unit(0.01,"mm"),legend.key.size = unit(3, "mm"))

PlotGenelogFacrossDay2

grid.arrange(lgddaysDEG_grob)
#facet5.5_ld<-grid.arrange(blank,lgddaysDEG_UpDown_grob, ncol=2,widths =c(7,1.5))


png(paste0(ofig,"Figure2_DEG_Selenop.png"),units = "in", width=10, height= 12, res = 600, family = "Arial")
ggpubr::ggarrange(facet5  , facet6_ld,facet5.5_ld,PlotGenelogFacrossDay2, ncol=2,nrow = 2, labels=c("A","B","C","D"))
dev.off()
