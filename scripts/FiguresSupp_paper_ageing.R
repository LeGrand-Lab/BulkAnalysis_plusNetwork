library(DESeq2)
library(rtracklayer)
library(dplyr)
library(tidyverse)
library(ggrepel)
library(WGCNA)
library(biomaRt)
library("ggsci")
library(cowplot)
library(apeglm)
library(ggplot2)
library(RColorBrewer)
library(ggthemes)
library(MASS)
library(pheatmap)
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
library(ggraph)
library(igraph)
library(viridis)
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

tableGSEA<-"GSEA/TableGSEA/"
fullGSEAconcatfile<-paste0(odir,tableGSEA,"GSEA_table_static_full.csv" )
fullGSEAconcat<-read.csv(fullGSEAconcatfile,sep = " " )
GSEAsigni<-fullGSEAconcat %>% filter(padj<=0.05) %>% mutate(sens=ifelse(NES >0,"Up","Down"))
infoGSEA<-readRDS(paste0(odir,tableGSEA,"TableGSEAsigniDynamics.rds"))

genes_df <- read.table("data/genesinfo.csv", sep="\t", header=T)
daysv<-list()
daysv<-lapply( unique(fullDEsta$type), function(t) unique(fullDEsta[fullDEsta$type == t,]$day) )
names(daysv)<-unique(fullDEsta$type)
Typecellv<-unique(fullDEsta$type)
days<-unique(fullDEsta$day)


orderTypecell=c("ECs","FAPs","MuSCs","Neutrophils","Inflammatory-Mac","Resolving-Mac")
colorsType=c("#10b387ff","#3d85c6ff","#b171f1ff","#f0e442ff","#ff9900ff","#cc0000ff")
names(colorsType) = orderTypecell
myannot = data.frame(CellType=factor(metadata$type, levels = orderTypecell),
                     Time = factor(metadata$time,levels= str_sort(unique(metadata$time))),
                     Age = factor(metadata$age, levels = str_sort(unique(metadata$age), decreasing =T)))
rownames(myannot) = rownames(metadata)

colorsTime = c(brewer.pal(9,"BuPu")[3],brewer.pal(9,"BuPu")[5],brewer.pal(9,"BuPu")[7],brewer.pal(9,"BuPu")[9])
names(colorsTime) = sort(unique(myannot$Time))
colorsAge = c("#ffc35d","#019190")
names(colorsAge) = c("Young","Old")
mycolors = list("CellType"=colorsType,"Time" = colorsTime, "Age" =colorsAge )


fullDEsta <- fullDEsta %>% mutate (type = factor( fullDEsta$type , levels =orderTypecell ))

# set aesthetics data:
sortie_full <- fullDEsta %>% mutate(DEsignificant=case_when(
  padj <= 0.05 & log2FoldChange >= 1.2 ~ "Signi UP" ,
  padj <= 0.05 & log2FoldChange <= -1.2 ~ "Signi DOWN",
  TRUE  ~ "Not Signi"
))

g <-ggplot(sortie_full, aes(x=log2FoldChange, y = -log10(padj+1e-20),color=DEsignificant)) +
  geom_point(aes(color=DEsignificant),size=.3) +
  scale_color_manual(values=c("lightgrey",brewer.pal(10, "RdBu")[9],brewer.pal(10, "RdBu")[2] )) +
  geom_vline(xintercept = c(1.2,-1.2), data=,color= "black",
             linetype="dashed", size=.2) +
  facet_grid_blank(vars(type),vars(day), drop = FALSE) + 
  theme_calc() +
  theme(panel.grid.major=element_blank()) +
  
  geom_text_repel(
    data= subset(sortie_full, padj < 0.0005 & 
                   abs(log2FoldChange) > 1.2),
    aes(label=symbol, fill=DEsignificant),
    size=2,
    segment.size = .1,
    force=2, force_pull = 2,
    max.overlaps=15
  ) +
  labs(title="Old vs Young across day&type", subtitle = "Signi = Differential gene padj < 0.05, -1.2<log2FoldChange>1.2",
       caption="vertical lines: ABS(log2FoldChange)=1.2 
          labels only for genes padj < 0.0005")+
  theme(axis.title = element_text(size=5, face = "bold"),axis.text=element_text(size=3),legend.title=element_text(size=6, face = "bold"),legend.text=element_text(size=4),title=element_text(size=7, face = "bold"),
        legend.spacing= unit(0.1,"points"),legend.spacing.y=unit(0.01,"points"),strip.text.x = element_text(size = 6, face = "bold", color="white"),strip.text.y = element_text(size=6, face = "bold", color="white"),legend.key.size = unit(3, "mm"),panel.spacing = unit(0.5,"mm"),panel.grid =element_line(color="white"))

#Add color in strip background
g2 <- ggplot_gtable(ggplot_build(g))
stripRowName <- which(grepl('strip-', g2$layout$name))
k <- 1
fills <- c(colorsTime,colorsType)
for (i in stripRowName) {
  j <- which(grepl('rect', g2$grobs[[i]]$grobs[[1]]$childrenOrder))
  g2$grobs[[i]]$grobs[[1]]$children[[j]]$gp$fill <- fills[k]
  k <- k+1
}
png(paste0(odir,"PlotsDEG/volcano_static.png"),units = "in", width=10, height= 5.5, res = 600, family = "Arial")
grid::grid.draw(g2)
dev.off()
tiff(paste0(odir,"PlotsDEG/volcano_static.tiff"),units = "in", width=10, height= 4, res = 600, family = "Arial")
grid::grid.draw(g2)
dev.off()
pdf(paste0(odir,"PlotsDEG/volcano_static.pdf"), width=14, height=10)
grid::grid.draw(g2)
dev.off()


# Hierarchical pathway enrichment overview
explorer<- function(hierachicalPathway2,node2,path,listSommetparcouru,data_edged, level){
  level=level+1
  listSommetparcouru<-c(listSommetparcouru,node2)
  newpath=paste0(path,".",node2)
  NamePath<-str_replace_all(node2,pattern = "_", replacement = " ")
  #size=GSEAtempo %>% filter(pathway==paste0("REACTOME_",node2)) %>% select(size) %>% unlist() %>% as.character() %>% as.numeric()
  size=ifelse(node2== "root1", 1,GSEAsigni %>% filter(pathway==paste0("REACTOME_",node2)) %>% dplyr::select(size) %>% unlist() %>% as.character() %>% as.numeric() %>% sum(na.rm = T))
  data_edged<-rbind(data_edged,c(path,newpath,level,size,NamePath))
  path=newpath
  for ( node2 in hierachicalPathway2[hierachicalPathway2$Node1 == node2,2]){
    if (!(node2 %in% listSommetparcouru)){
      data_edged<-explorer(hierachicalPathway2,node2,path,listSommetparcouru,data_edged,level)
    }
  }
  return(data_edged)
}  

ConvertPathwaysName<-read.csv("data/ReactomePathways.txt",sep = "\t",header = F)
ConvertPathwaysName <- ConvertPathwaysName %>% filter(str_detect(V3,"Mus musculus"))

hierachicalPathway<-read.csv("data/ReactomePathwaysRelation.txt",sep = "\t",header = F)
hierachicalPathway<-hierachicalPathway %>% filter(str_detect(V1,"R-MMU")) %>% mutate("Node1"= str_replace_all(str_to_upper(ConvertPathwaysName[na.omit(match(hierachicalPathway$V1,ConvertPathwaysName$V1)),]$V2), "[\\ /\\-\\:]",'_'))  %>% mutate("Node1"= str_replace_all(Node1,"[,()]", ""))
hierachicalPathway<-hierachicalPathway %>% filter(str_detect(V2,"R-MMU")) %>% mutate("Node2"= str_replace_all(str_to_upper(ConvertPathwaysName[na.omit(match(hierachicalPathway$V2,ConvertPathwaysName$V1)),]$V2),"[\\ /\\-\\:]",'_')) %>% mutate("Node2"= str_replace_all(Node2,"[,()]", ""))

topHierarchie<-setdiff(hierachicalPathway$Node1,hierachicalPathway$Node2)
hierachicalPathway2<- rbind(cbind(hierachicalPathway$Node1,hierachicalPathway$Node2),cbind(rep("root1",by=length(topHierarchie)),topHierarchie))

colnames(hierachicalPathway2)<-c("Node1","Node2")
hierachicalPathway2<-as.data.frame(hierachicalPathway2)

listSommetparcouru<-c()
data_edged=data.frame()
data_edged<-explorer(hierachicalPathway2,"root1","root1",listSommetparcouru,data_edged,0)
colnames(data_edged)<-c("Node1","Node2","level","size","NamePath")
data_edged<- data_edged %>% mutate(level=factor(as.numeric(level),levels = 1:13)) 

edges <- data_edged %>% dplyr::select(Node1,Node2)
minNES<-min(infoGSEA$NES)
maxNES<-max(infoGSEA$NES)
listCircles<-list()
for( t in Typecellv){
  for(d in daysv[[t]]){
    print(paste0(t,'_',d))
    infoGSEA2<- infoGSEA %>% dplyr::filter(day==d,type==t)
    
    vertices <- data_edged%>% mutate(showlabel=ifelse(level == 2 , NamePath,NA)) %>% dplyr::select (Node2,size,NamePath,level,showlabel) %>% rbind(c("root1",7,"root",1,NA))
    vertices[vertices$NamePath == "REPRODUCTION",]$showlabel<-""
    vertices <- vertices %>% mutate(NbCondSigni=infoGSEA2[match(paste0(" ", vertices$NamePath), infoGSEA2$NamePath),]$NbCondSigni) %>% mutate(NES=infoGSEA2[match(paste0(" ", vertices$NamePath), infoGSEA2$NamePath),]$NES)  %>% mutate(NbCondSigni=ifelse(is.na(NbCondSigni), 0, NbCondSigni)) %>% mutate(NES=ifelse(is.na(NES), 0, NES))
    verticesHeader<-vertices %>% dplyr::filter(as.numeric(level)>=2) %>% dplyr::select(NamePath,NbCondSigni) %>% mutate(name=str_replace_all(NamePath," ","_")) %>% mutate(name=str_replace_all(name,"\\+","\\\\\\+"))
    verticesHeader<- verticesHeader %>% mutate(NbTotPathEnrich=lapply(verticesHeader$name , function(n) nrow(vertices[grep(paste0("\\.",n,"\\."),vertices$Node2),] %>% filter(NbCondSigni>0))) %>% unlist()) %>% mutate(NbTotPathEnrich=NbTotPathEnrich+ifelse(NbCondSigni>0,1,0)) %>% filter(NbTotPathEnrich>0)
    
    vertices <- vertices %>% mutate(NbTotPathEnrich=verticesHeader[match(vertices$NamePath,verticesHeader$NamePath),]$NbTotPathEnrich)  %>% mutate(level2=ifelse(is.na(NbTotPathEnrich), NA,level)) %>% mutate(NbTotPathEnrich=ifelse(is.na(NbTotPathEnrich), 0,NbTotPathEnrich))
    vertices[vertices$NamePath =="root1",]$level2<-1
    vertices[vertices$NamePath =="root1",]$NbTotPathEnrich<-1
    vertices[vertices$NamePath =="root",]$level2<-1
    vertices[vertices$NamePath =="root",]$NbTotPathEnrich<-1
    vertices<-vertices %>% filter(NbTotPathEnrich>0)
    edge2<-edges %>% filter(Node2 %in% vertices$Node2)
    
    mygraph <- graph_from_data_frame( edge2, vertices=vertices )
    

    circles<-ggraph(mygraph, layout = 'circlepack', weight=as.numeric(NbTotPathEnrich)) + 
      geom_node_circle(aes(fill= NES ,group = level2)) +
      geom_node_label(aes(label=showlabel,size=NbTotPathEnrich),repel = T) +
      theme_void() + ggtitle(paste0(t, " ", d)) +
      theme(legend.position=c(fill="bottom"),legend.title = element_text(size=3,face="bold"),legend.text=element_text(size=2), legend.key.size = unit(0.1,'cm'), plot.title=element_text(face="bold",size=10,hjust=0.5)) + 
      scale_fill_gradient2(name= "Number of celltype_day \n with pathway enriched",limits=c(minNES,maxNES) ,low = muted("blue"),mid="white",high=muted("red"), midpoint = 0) + scale_size_continuous(name="Number pathway enriched \n in this event hierarchy",range= c(0.5,2))
    
    #print(circles)
    listCircles[[paste0(t,'_',d)]]<-circles
    
    #png(file=paste0(ofig,"HirarchieALLpathwaysNEs_",t,"_",d,".png"),units = "in", width=8, height= 8, res = 800, family = "Arial")
    #print(circles)
    #dev.off()
       

    
  }
}
blank<-grob()
png(paste0(ofig,"HirarchieALLpathwaysNEs_UP_DOWN1.png"),units = "in", width=10, height= 15, res = 600, family = "Arial")
ggpubr::ggarrange(ggpubr::ggarrange(listCircles[["ECs_D0"]],listCircles[["ECs_D2"]],listCircles[["ECs_D4"]],listCircles[["ECs_D7"]], ncol = 4),
                  ggpubr::ggarrange(listCircles[["FAPs_D0"]],listCircles[["FAPs_D2"]],listCircles[["FAPs_D4"]],listCircles[["FAPs_D7"]], ncol = 4),
                  ggpubr::ggarrange(listCircles[["MuSCs_D0"]],listCircles[["MuSCs_D2"]],listCircles[["MuSCs_D4"]],listCircles[["MuSCs_D7"]], ncol = 4),
                  ggpubr::ggarrange(blank,listCircles[["Neutrophils_D2"]],blank,blank, ncol = 4),
                  ggpubr::ggarrange(blank,listCircles[["Inflammatory-Mac_D2"]],listCircles[["Inflammatory-Mac_D4"]],blank, ncol = 4),
                  ggpubr::ggarrange(blank,listCircles[["Resolving-Mac_D2"]],listCircles[["Resolving-Mac_D4"]],listCircles[["Resolving-Mac_D7"]], ncol = 4),nrow=6)
dev.off()


tiff(paste0(ofig,"HirarchieALLpathwaysNEs_UP_DOWN.tiff"),units = "in", width=10, height= 15, res = 600, family = "Arial")
ggpubr::ggarrange(ggpubr::ggarrange(listCircles[["ECs_D0"]],listCircles[["ECs_D2"]],listCircles[["ECs_D4"]],listCircles[["ECs_D7"]], ncol = 4),
                  ggpubr::ggarrange(listCircles[["FAPs_D0"]],listCircles[["FAPs_D2"]],listCircles[["FAPs_D4"]],listCircles[["FAPs_D7"]], ncol = 4),
                  ggpubr::ggarrange(listCircles[["MuSCs_D0"]],listCircles[["MuSCs_D2"]],listCircles[["MuSCs_D4"]],listCircles[["MuSCs_D7"]], ncol = 4),
                  ggpubr::ggarrange(blank,listCircles[["Neutrophils_D2"]],blank,blank, ncol = 4),
                  ggpubr::ggarrange(blank,listCircles[["Inflammatory-Mac_D2"]],listCircles[["Inflammatory-Mac_D4"]],blank, ncol = 4),
                  ggpubr::ggarrange(blank,listCircles[["Resolving-Mac_D2"]],listCircles[["Resolving-Mac_D4"]],listCircles[["Resolving-Mac_D7"]], ncol = 4),nrow=6)
dev.off()

pdf(paste0(ofig,"HirarchieALLpathwaysNEs_UP_DOWN.pdf"),units = "in", width=10, height= 14, res = 600, family = "Arial")
ggpubr::ggarrange(ggpubr::ggarrange(listCircles[["ECs_D0"]],listCircles[["ECs_D2"]],listCircles[["ECs_D4"]],listCircles[["ECs_D7"]], ncol = 4),
                  ggpubr::ggarrange(listCircles[["FAPs_D0"]],listCircles[["FAPs_D2"]],listCircles[["FAPs_D4"]],listCircles[["FAPs_D7"]], ncol = 4),
                  ggpubr::ggarrange(listCircles[["MuSCs_D0"]],listCircles[["MuSCs_D2"]],listCircles[["MuSCs_D4"]],listCircles[["MuSCs_D7"]], ncol = 4),
                  ggpubr::ggarrange(blank,listCircles[["Neutrophils_D2"]],blank,blank, ncol = 4),
                  ggpubr::ggarrange(blank,listCircles[["Inflammatory-Mac_D2"]],listCircles[["Inflammatory-Mac_D4"]],blank, ncol = 4),
                  ggpubr::ggarrange(blank,listCircles[["Resolving-Mac_D2"]],listCircles[["Resolving-Mac_D4"]],listCircles[["Resolving-Mac_D7"]], ncol = 4),nrow=6)
dev.off()

#tableSankeyplot function to extrat table nodes sources and target with value , for a sankey plot
#plot gathered gene with the same pattern (for one type cell, witch days are DEG and in witch sens) and links days in function of sens of log2foldchange
#so if gene is UP, it is more exessed in OLD condition than Young
#signiDEgene = table with all DEG padj<0.05
#vectTypeCell = type cell selected, vector possible
#minCombi = founds genes with at least 1 or 2 days DEG
#NODIFF = Nodes inclus days UNDIFF

tableSankeyplot<- function(signiDEgene,vectTypeCell,minCombi,NODIFF){
  sourceNode<-c()
  targetNode<-c()
  valueEdge<-c()
  NbUpDown<-c()
  Nbconnection<-c()
  CombiDayDEG<-c()
  
  for ( typeCell in vectTypeCell ) {
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
  links <- data.frame(
    from=links1$from,
    to=links1$to,
    substance= as.factor(links1 %>% dplyr::select(group) %>% unlist()),
    #substance= as.factor(links1$CombiDayDEG),
    quantity=links1$quantity
  )
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
  nodes<-nodes  %>% arrange(loops,sens,day)
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
  my_title <- outputprefix
  attr(my_title, "gp") <- grid::gpar(fontsize=12, fontface="bold", col="black")
  
  # node style
  ns <- list(type="arrow",gp=gpar(fill=dblue, col="white", lwd=3),
             length=0.5,
             label_gp=gpar(col="black", fontsize=6,fontface="bold"),
             mag_pos="label", mag_fmt="%.0f", mag_gp=gpar(fontsize=6,fontface="bold",col=dblue))
  
  #pdf(paste0(analyse,"test_group_",group,"_",titleNODIFF,'_',sqrt,'_',outputprefix,"_sankeyColorUPDOWN.pdf"), width=10, height=7) # Set up PDF device
  PantaRhei::sankey(nodes, links, palettes,
                    max_width=0.15, rmin=0.5,
                    node_style=ns,
                    page_margin=c(0, 0.05, 0, 0.05),
                    title=my_title, legend=gpar(fontsize=4, col="blue", ncols=1))
  #dev.off()
  
}

#Vectcorlor in function group chosen
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


dataFrameTotalUPDOWN<-readRDS(paste0(pathcombinaisontableUPDOWN,'.rds'))
orderTypecell=c("ECs","FAPs","MuSCs","Neutrophils","Inflammatory-Mac","Resolving-Mac")
colorsType=c("#10b387ff","#3d85c6ff","#b171f1ff","#f0e442ff","#ff9900ff","#cc0000ff")
names(colorsType)=orderTypecell
orderTypecell=c("ECs","FAPs","MuSCs","Inflammatory-Mac","Resolving-Mac")
for(typecell in orderTypecell){
tableSankeyplotFAPs<-tableSankeyplot(signiDEgene,typecell,1,F)
png(file=paste0(ofig,"SigniGeneFlow_",typecell,".png"),units = "in", width=6, height= 4, res = 800, family = "Arial")
sankeyPantaRhei(tableSankeyplotFAPs,colorsType[[typecell]],"NbUpDown",colorNbUpDown,typecell,T,F)
dev.off()
}

#Heatmap
library(msigdbr)
thegmt <- msigdbr(species = "Mus musculus", 
                  category = 'C2', 
                  subcategory=c('CP:REACTOME'))
msigdbr_list = split(x = thegmt$gene_symbol, f = thegmt$gs_name)

for( t in Typecellv){
  for(d in daysv[[t]]){
    print(paste0(t,'_',d))
    infoGSEA2<- infoGSEA %>% dplyr::filter(day==d,type==t) %>% arrange(padj) %>% dplyr::select(NbCondSigni,pathway,NES,DayTypeSens,leadingEdge) %>% head()
    genes<-infoGSEA2 %>% dplyr::select(leadingEdge)
    print(genes)
    infogenes<-lapply(infoGSEA2$pathway ,function(p) signiDEgene %>% filter(day==d,type==t,symbol %in% c(msigdbr_list[[p]]) ))
  }
}    
    


