# Explore result GSEA statics,
# 1. Observing the patterns (Day, Cell type, Sens) of differential enriched pathways
# 2. Examined pathway with heatmap, plot of expression quantity 
# 3. Organised pathway with hierarchies of Reactome
# 
# --
# PaulineM
library(tidyverse)
library(ComplexHeatmap)
library(msigdbr)
library(cowplot)
library(openxlsx)
library(arcdiagram)
library(devtools)
library(igraph)
library(DESeq2)
library(ggrepel) # for labels to point
library(rmarkdown)
library(ggnewscale)
library(reshape2)
library(circlize)
library(oce)
library(ggraph)
library(arcdiagram)
library(plyr)
library(circlepackeR) 
library(igraph)
library(data.tree)
library(spaceNtime)
library(RColorBrewer)
library(patchwork)
library(facetious)
##################################
#### Path files to load or created
##################################

#Directories
setwd("~/BulkAnalysis_plusNetwork/")
dir<-"/home/bioinfo/BulkAnalysis_plusNetwork/"
odir <- paste0(dir,"/exam_INTER_conditions/static/")

plots<-"GSEA/PlotsGSEA/"
tableGSEA<-"GSEA/TableGSEA/"
NormData <- "data/CountNormalised/"
HierarchieData= "GSEA/HierarchieData/"

#Loads

gseaoutfull2padjfile = paste0(odir,tableGSEA,"GSEA_table_static_full.rds" )
fullGSEAconcatfile<-paste0(odir,tableGSEA,"GSEA_table_static_full.csv" )
gseaout_filtered_full2padjfile = paste0(odir,tableGSEA,"GSEA_table_static_softfilter.rds" )
countnormalisedfile=paste0(NormData,"CountNormalised.txt")
meancountnormalisedfile=paste0(NormData,"MeanCountNormalised.txt")
DEGtableTotalUPDOWNfile<-paste0(odir,"TableDEG/TableDynamicUpDownDEG.rds")
DEGtablefile<-paste0(odir,"TableDEG/DEG_table_static_full.rds")


#Created
#Completed result GSEA with pattern associated to pathway 
TableGSEADynamics<-"TableGSEAsigniDynamics"
plotpaways<-paste0(plots,"Pathways/")
plotpathwaysHierarchie<-paste0(plots,"PathwaysHierarchie/")
reportsuppPathway<-"reports/reportPathwaysEnrichedStatics/"
reportsuppHierarchie<-"reports/reportPathwaysHierarchie/"
listHierarchieREACTOME<-paste0(odir,HierarchieData,"ListReactomeHiearachieUseFunction.txt")
tablePathwayHierarchieFull<-paste0(odir,HierarchieData,"TableReactomeHiearachieUseFunction.txt")
tablePathwayHierarchieTotalPathwayEnrichie<-paste0(odir,HierarchieData,"TableReactomeHiearachieTotalPathwayEnrichie.txt")


##################################
#### Load files
##################################
#rerun code to created y=table and plot

fullGSEA = readRDS(gseaoutfull2padjfile)
fullGSEAconcat<-read.csv(fullGSEAconcatfile,sep = " " )
filteredGSEA =readRDS(gseaout_filtered_full2padjfile) 
countnormalised=read.table(countnormalisedfile,header = T)
meancountnormalised=read.table(meancountnormalisedfile,header = T)
DEGtableTotalUPDOWN<-readRDS(DEGtableTotalUPDOWNfile)
DEGtable<-readRDS(DEGtablefile)

GSEAsigni<-fullGSEAconcat %>% filter(padj<=0.05) %>% mutate(sens=ifelse(NES >0,"Up","Down"))

#Extract patway and gene associated , same use for fgsea analysis
thegmt <- msigdbr(species = "Mus musculus", 
                  category = 'C2', 
                  subcategory=c('CP:REACTOME'))
msigdbr_list = split(x = thegmt$gene_symbol, f = thegmt$gs_name)

genes_df <- read.table("data/genesinfo.csv", sep="\t", header=T)

daysv<-list()
daysv<-lapply( unique(fullGSEAconcat$type), function(t) unique(fullGSEAconcat[fullGSEAconcat$type == t,]$day) )
names(daysv)<-unique(fullGSEAconcat$type)
Typecellv<-unique(fullGSEAconcat$type)

days<-unique(fullGSEAconcat$day)

####################
####Function colors
####################

orderTypecell=c("ECs","FAPs","MuSCs","Neutrophils","Inflammatory-Mac","Resolving-Mac")
colorsType=c("#10b387ff","#3d85c6ff","#b171f1ff","#f0e442ff","#ff9900ff","#cc0000ff")
hsl_color<-data.frame("ECs"=164,"FAPs"=208,"MuSCs"=270,"Neutrophils"=56,"Inflammatory-Mac"=36,"Resolving-Mac"=0)
names(colorsType)=orderTypecell
colorsTime = c(brewer.pal(9,"BuPu")[3],brewer.pal(9,"BuPu")[5],brewer.pal(9,"BuPu")[7],brewer.pal(9,"BuPu")[9])
names(colorsTime)<-days
colorsAge = c("#ffc35d","#019190")
names(colorsAge) = c("Young","Old")
orderDayAge=lapply(days, function(d) c(paste(d,"Young"),paste(d,"Old"))) %>% unlist()

################################################# 
# Extract stats significant pathways enrichied
#################################################
NumberPathwaysSigni = array(NA, dim=c(length(days),length(Typecellv)))
rownames(NumberPathwaysSigni) =  days
colnames(NumberPathwaysSigni) = Typecellv

NumberPathwaysUniqueOnOneDay<-c()
NumberPathwaysUniqueOnOneTypeCell<-c()
for (d in days){
  NumberPathwaysUniqueOnOneDay<-c(NumberPathwaysUniqueOnOneDay,GSEAsigni %>% filter(day==d) %>% dplyr::select(pathway) %>% unique() %>% unlist() %>% length())
  for (t in Typecellv){
    NumberPathwaysSigni[d, t] <- dim(GSEAsigni %>% filter(day==d & type==t))[1]
    if (d == days[1]){
      NumberPathwaysUniqueOnOneTypeCell<-c(NumberPathwaysUniqueOnOneTypeCell, GSEAsigni %>% filter(type==t) %>% dplyr::select(pathway) %>% unique() %>% unlist() %>% length())
    }
  }
}

names(NumberPathwaysUniqueOnOneDay)<-days
names(NumberPathwaysUniqueOnOneTypeCell)<-Typecellv
NumberPathwaysSigni
#Plot 

Violindata<-data.frame(ViolinDT=paste0(GSEAsigni$day,"_",GSEAsigni$type),
                       ViolinDay=factor(GSEAsigni$day,levels= unique(GSEAsigni$day)),
                       CellType=factor(GSEAsigni$type, levels= orderTypecell),
                       ViolinValue=GSEAsigni$NES,
                       nbPathwaySignificant = lapply(1:length(GSEAsigni$NES), function(i) NumberPathwaysSigni[GSEAsigni[i,]$day,GSEAsigni[i,]$type] ) %>% unlist(),
                       a=rep("a",length(GSEAsigni$NES)))


PlotDEGonlyOneDayVioline<-ggplot(Violindata, aes(x=a ,y=ViolinValue, fill=nbPathwaySignificant)) +# fill=name allow to automatically dedicate a color for each group
  geom_jitter( shape=16, size=1, position=position_jitter(0.4), aes(color=nbPathwaySignificant,alpha=0.1), show.legend=F)+
  geom_violin(scale="count",aes(alpha=0.01),show.legend = c("alpha"=F))+
  facet_grid_blank(vars(CellType),vars(ViolinDay), drop = FALSE) +
  scale_fill_gradient(low = "#E7E1EF",high="#67001F",labels = scales::label_comma())+
  scale_color_gradient(low = "#E7E1EF",high="#67001F",labels = scales::label_comma())+
  coord_flip()+
  scale_x_discrete(limits=rev)  +   
  geom_hline(yintercept = 0, color="black", linetype="dashed",lwd=0.5)+
  ylab("log2FoldChange")+xlab("Type Cell") +
  theme(axis.title = element_text(size=5, face = "bold"),axis.text.x=element_text(size=3,colour = "black", face="bold"),axis.text.y=element_text(colour = "white"),legend.title=element_text( size=6, face = "bold"),legend.text=element_text(size=4),title=element_text(size=7, face = "bold"),
        legend.spacing= unit(0.1,"points"),legend.position = "bottom",legend.spacing.y=unit(0.01,"mm"),legend.key.size = unit(3, "mm"),panel.grid =element_line(color="white"),strip.text.x = element_text(color = "white", face= "bold",size=7),strip.text.y = element_text(color = "black", face= "bold",size=7),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.background = element_blank())+labs(fill = "Count DEG padj<0.05")+
  ggtitle("Distribution of significatif differentially expressed genes Old vs Young \n according to log2Foldchange")
#Add color in background grid
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
png(paste0(odir,plots,"PlotPathwaysigniDayTypeVioline.png"),units = "in", width=5, height= 7, res = 300, family = "Arial")
grid::grid.draw(g5)
dev.off()

#### Explore Dynamique pathways Enrichied padj<0.05

## table of number of pathways found in 1, 2 , ... conditions
RepartitionPathwaySummary<-table(table(GSEAsigni$pathway))
## dataframe : Var1 = pathways Freq = Number condition which the pathway is found
RepartitionPathway<-as.data.frame(table(GSEAsigni$pathway)) 

## NewTable with NbCondSigni, CondSigni: Day1_Type1_Sens1.Day_Type2_Sens2 ..., DayTypeSens,NbPathWithDynamique,color, links to report paway
DynamicPathwayAcrossCond <- function(RepartitionPathway,GSEAsigni,n){
  PathwaySigniCond<-RepartitionPathway %>% filter(Freq==n)
  PathwaySigniCondtable<- lapply(PathwaySigniCond$Var1 , function(x) dplyr::filter(GSEAsigni, pathway == x))
  DynamicCond <-lapply(PathwaySigniCondtable, function(x) paste0(x$day,'_',x$type,'_',x$sens) %>% sort())
  names(DynamicCond)<-lapply(PathwaySigniCondtable, function(x) unique(x$pathway) )
  DynamicCondconcat<-c()
  for(Dym in DynamicCond){DynamicCondconcat=c(DynamicCondconcat,str_c(Dym,collapse = "."))}
  names(DynamicCondconcat)<-lapply(PathwaySigniCondtable, function(x) unique(x$pathway))
  Dynamiquetable<-ldply(lapply(PathwaySigniCondtable, function(x) as.data.frame(cbind(rep(n,n),x,paste0(x$day,'_',x$type,'_',x$sens),rep(DynamicCondconcat[[unique(x$pathway)]],n)))),rbind)
  colnames(Dynamiquetable)<-c("NbCondSigni",colnames(PathwaySigniCondtable[[1]]),"DayTypeSens","CondSigni")
  
  Dynamique<-lapply(unique(Dynamiquetable$pathway), function(x) unique(Dynamiquetable %>% filter(pathway==x) %>% dplyr::select(CondSigni) %>% unlist())) %>% unlist()
  names(Dynamique)<-unique(Dynamiquetable$pathway)
  tableDynamique<-table(Dynamique)
  Dynamiquetable <- Dynamiquetable %>% mutate(NbPathWithDynamique=as.numeric(as.character(tableDynamique[match(Dynamiquetable$CondSigni,names(tableDynamique))]))) %>% mutate(color=lapply(Dynamiquetable$type, function(t) colorsType[[t]]) %>% unlist()) %>% arrange(NbPathWithDynamique,decreasing =T)
  Dynamiquetable <- Dynamiquetable %>% mutate( DayType = paste0(day,'_',type))
  Dynamiquetable <-  Dynamiquetable %>% mutate( link=paste0("<a href='../",reportsuppPathway,"PathwayEnrichedIn",n,"DayTypeCell_",pathway,"_meanCountNormalizedGene.html","'>FigureSupp</a>"))
  Dynamiquetable <-  Dynamiquetable %>% mutate( NamePathway= str_replace(str_replace_all(pathway,"_"," "), "REACTOME",""))
  return(Dynamiquetable)
}

Dynamic8table=DynamicPathwayAcrossCond(RepartitionPathway,GSEAsigni,8)
Dynamic7table=DynamicPathwayAcrossCond(RepartitionPathway,GSEAsigni,7)
Dynamic6table=DynamicPathwayAcrossCond(RepartitionPathway,GSEAsigni,6)
Dynamic5table=DynamicPathwayAcrossCond(RepartitionPathway,GSEAsigni,5)
Dynamic4table=DynamicPathwayAcrossCond(RepartitionPathway,GSEAsigni,4)
Dynamic3table=DynamicPathwayAcrossCond(RepartitionPathway,GSEAsigni,3)
Dynamic2table=DynamicPathwayAcrossCond(RepartitionPathway,GSEAsigni,2)
Dynamic1table=DynamicPathwayAcrossCond(RepartitionPathway,GSEAsigni,1)

#Plot of Pathway in function NES
GetPlot<-function(Dynamic6table,GSEAsigni, nbCond , vertical){
  Dynamic6table <-Dynamic6table %>% mutate( CondSigni = as.factor(CondSigni)) %>% mutate(CondSigni = fct_reorder(CondSigni, as.numeric(NbPathWithDynamique), .desc = T) )
  Dynamic6table<-Dynamic6table %>% mutate(type = factor(Dynamic6table$type,orderTypecell ))
  col<- Dynamic6table %>% arrange(type) %>% dplyr::select(color) %>% unlist() %>% as.character() %>% unique()
  
  typebytime<-ggplot(Dynamic6table, aes(x=day,y=as.numeric(NES),fill=type,color=type,group=paste0(pathway,'_',type))) +
    geom_point(aes(size=as.numeric(NbPathWithDynamique),alpha=0.7),show.legend = c(size=T,alpha=F,col=F,fill=F))+
    scale_size_binned(name="Number path",limits=c(1,MaxNbPaathWithSameDyn))+
    ylab("NES")+
    geom_line()+
    geom_hline(yintercept = 0, color="black", linetype="dashed",lwd=0.5,slope=F)+
    scale_color_manual(name="Cell type",values=col)+
    ggtitle( paste0("Pathways with ",nbCond," conditions enriched")) +
    facet_wrap(~CondSigni+NbPathWithDynamique ,dir = "v", ncol=vertical,labeller = labeller(NbPathWithDynamique = label_facet(Dynamic6table$NbPathWithDynamique,"Number pathways: ")))+
    theme_linedraw()+ theme(strip.text  = element_text(size = 11 , face = "bold"),strip.background = element_rect(fill="grey"),panel.grid.major.x = element_line(color = "gray60", size = 0.8),axis.title= element_text(face="bold"),axis.text = element_text(face="bold",size=14),legend.position="bottom")
  print(typebytime)
  save_plot(paste0(odir, plots, 'Pathway_enriched_in__',nbCond ,'conditions_Plot2.svg'),typebytime, base_height = 10,
            base_width = 14)
}



MaxNbPaathWithSameDyn<-max(Dynamic8table$NbPathWithDynamique,Dynamic7table$NbPathWithDynamique,Dynamic6table$NbPathWithDynamique,Dynamic5table$NbPathWithDynamique,Dynamic4table$NbPathWithDynamique,Dynamic3table$NbPathWithDynamique,Dynamic2table$NbPathWithDynamique,Dynamic1table$NbPathWithDynamique )

GetPlot(Dynamic8table,GSEAsigni, 8 , 1)
GetPlot(Dynamic7table,GSEAsigni, 7 , 1)
GetPlot(Dynamic6table,GSEAsigni, 6 , 1)
GetPlot(Dynamic5table,GSEAsigni, 5 , 2)
GetPlot(Dynamic4table,fullGSEAconcat, 4 , 2)
GetPlot(Dynamic3table,GSEAsigni, 3 , NA)
GetPlot(Dynamic2table,GSEAsigni, 2 , NA)
GetPlot(Dynamic1table,GSEAsigni, 1 , NA)

##############
# Explore data Expression associated to dynamic tendance found in enriched pathway
##############
minCountNormalized<-log2(round(min(as.numeric(DEGtableTotalUPDOWN$meanCountNormalizedYoung,DEGtableTotalUPDOWN$meanCountNormalizedOld, na.rm=T )),digits = 0)+1)
maxCountNormalized<-log2(round(max(as.numeric(DEGtableTotalUPDOWN$meanCountNormalizedYoung,DEGtableTotalUPDOWN$meanCountNormalizedOld, na.rm=T )),digits = 0)+1)
  
minlogF <-round(min(DEGtableTotalUPDOWN$log2FoldChange, na.rm=T ),digits = 0)
maxlogF <- round(max(DEGtableTotalUPDOWN$log2FoldChange, na.rm=T ),digits = 0)

DEGtableTotalUPDOWN<- DEGtableTotalUPDOWN %>% unite( 'GeneDayType' , c(symbol,day,type) ,sep = '_',remove = F)


matricePathHeatmap<-function(path, tablePlottest){
  genes_list<-lapply(unique(tablePlottest$DayType), function(x) tablePlottest %>% filter(age=="Young",DayType == x) %>% arrange(padj,log2F)  %>% filter(padj <= 0.05) %>% dplyr::select(symbol) %>% unlist() %>% unique() %>% as.character())
  genes<-genes_list %>% unlist() %>% as.character() %>% unique()
  genes<-na.exclude(genes)
  maxgene<-max(lapply(1:length(genes_list), function(i) length(genes_list[[i]])) %>% unlist())
  if (length(genes) > 40){
    while (length(genes) > 40 ){
      maxgene<-maxgene-1
      genes<-lapply(1:length(genes_list), function(i) genes_list[[i]][1:maxgene]) %>% unlist() %>% as.character() %>% unique()
      genes<-na.exclude(genes)
      
      }
  }
  matrix<-data.frame(symbol=genes )
  colNAmesMatrix<-c("symbol")
  typev<-c()
  for ( t in orderTypecell){
    for ( d in daysv[[t]]){
      colNAmesMatrix<-c(colNAmesMatrix,paste0(t,'_',d))
      typev<-c(typev,t)
      DEGtableTotalUPDOWN_tempo<-DEGtableTotalUPDOWN %>% filter(day==d,type==t)
      matrix <-matrix %>% mutate( colNAmesMatrix = as.numeric(DEGtableTotalUPDOWN_tempo[match(matrix$symbol,DEGtableTotalUPDOWN_tempo$symbol),]$log2FoldChange))
      colnames(matrix)<-colNAmesMatrix
      }
  }
  rownames(matrix)<-matrix$symbol
  matrix<-matrix[,-1]
  typev<-factor(typev,levels = orderTypecell)
  colgroup<-as.data.frame(typev)
  rownames(colgroup)<-colNAmesMatrix[2:length(colNAmesMatrix)]
  ha = HeatmapAnnotation(foo = anno_block(gp = gpar(fill =colorsType )))
  heatmap1<-ComplexHeatmap::Heatmap(as.matrix(matrix),cluster_columns = FALSE,cluster_rows = FALSE ,column_split = colgroup,top_annotation = ha,column_names_rot = 45,name=path,heatmap_legend_param = list( title="Log2FoldChange",
                                                                                                                                                                            direction = "horizontal"),use_raster = TRUE,col=colorRamp2(c(minlogF,0,maxlogF), c("blue", "white", "red")))
  heatmap2<-ComplexHeatmap::Heatmap(as.matrix(matrix),cluster_columns = FALSE,cluster_rows = FALSE ,column_split = colgroup,top_annotation = ha,column_names_rot = 45,name=path,heatmap_legend_param = list( title="Log2FoldChange",
                                                                                                                                                                                                            direction = "horizontal"),use_raster = TRUE,col=colorRamp2(c(min(matrix, na.rm=T),0,max(matrix, na.rm=T)), c("blue", "white", "red")))
  
  return(list(heatmap1,heatmap2))
  }


  
PlotGeneCountNorm <- function(Dynamic7table,msigdbr_list,numberCond){

  pathway_list<-unique(Dynamic7table$pathway)
  gene_list<-lapply(1:length(pathway_list), function(x) msigdbr_list[pathway_list[x]] %>% unlist() %>% unique())
  names(gene_list)<-pathway_list
  DynamicTableSmall<- Dynamic7table %>% dplyr::select(CondSigni,pathway,day,type,DayType,DayTypeSens,NES,color)
  for ( path in  pathway_list){
    if ( file.exists(paste0(dir,reportsuppPathway,"PathwayEnrichedIn",numberCond,"DayTypeCell_",path,"_meanCountNormalizedGene.html")) == F){
    tablePlot<-data_frame()
    numberGeneInPathway=length(unlist(gene_list[path]))
    DynamicPath<-DynamicTableSmall %>% filter(pathway==path) %>% mutate(type=factor(type,level=orderTypecell))
    NameDinamic<-unique(DynamicPath$CondSigni)
    tablePlotint<-cbind(symbol=unlist(gene_list[path]),do.call("rbind",replicate(numberGeneInPathway,DynamicPath,simplify = FALSE))%>%arrange(DayType))
    tablePlotint<-tablePlotint %>% mutate(GeneDayType=paste0(symbol,'_',day,"_",type))
    tablePlot<-rbind(tablePlot,tablePlotint)
    #DEGtableTotalUPDOWN<- DEGtableTotalUPDOWN %>% unite( 'GeneDayType' , c(symbol,day,type) ,sep = '_',remove = F)
    tablePlot<- tablePlot %>% mutate( padj = as.numeric(DEGtableTotalUPDOWN[match(tablePlot$GeneDayType,DEGtableTotalUPDOWN$GeneDayType),]$padj))
    tablePlot<- tablePlot %>% mutate( log2F = as.numeric(DEGtableTotalUPDOWN[match(tablePlot$GeneDayType,DEGtableTotalUPDOWN$GeneDayType),]$log2FoldChange))
    tablePlotY<- tablePlot %>% mutate( meanCountNormalized = as.numeric(DEGtableTotalUPDOWN[match(tablePlot$GeneDayType,DEGtableTotalUPDOWN$GeneDayType),]$meanCountNormalizedYoung)) %>% mutate( age = rep("Young",length(tablePlot$day)))
    tablePlotO<- tablePlot %>% mutate( meanCountNormalized = as.numeric(DEGtableTotalUPDOWN[match(tablePlot$GeneDayType,DEGtableTotalUPDOWN$GeneDayType),]$meanCountNormalizedOld)) %>% mutate( age = rep("Old",length(tablePlot$day)))
    tablePlot <- rbind(tablePlotY,tablePlotO) %>% filter(padj < 0.05)
    tablePlottest <- tablePlot %>% unite("DayAge",c(day, age),sep = ' ',remove = F) 
    tablePlottest <- tablePlottest %>% arrange(day,age,decreasing = c(F, T)) 
    tablePlottest <- tablePlottest %>% mutate (DayAge = factor( tablePlottest$DayAge , levels =orderDayAge ))
    
    geneDEGdt_list<-lapply(unique(DynamicPath %>% dplyr::select(DayType) %>% unlist()), function(x) tablePlottest %>% filter(DayType == x) %>% dplyr::select(symbol) %>% unlist() %>% as.character())
    names(geneDEGdt_list)<-unique(DynamicPath %>%  dplyr::select(DayType) %>% unlist())
    leadingEdge_list <- lapply(unique(DynamicPath %>% dplyr::select(DayType) %>% unlist()), function(x) Dynamic7table %>% filter(pathway == path,DayType == x) %>% dplyr::select(leadingEdge) %>% unlist() %>% str_split(pattern = " " ) %>% unlist() %>% str_replace_all(pattern = "\n", replacement = "") %>% keep( ~ any(geneDEGdt_list[[x]] %in% .x)))  
    leadingEdge<-data.frame(x=c(length(leadingEdge_list[[1]]), paste(leadingEdge_list[[1]],collapse = " ")))
    
    if(length(leadingEdge_list) >1){
      for( x in 2:length(leadingEdge_list)){leadingEdge<-cbind(leadingEdge,c(length(leadingEdge_list[[x]]), paste(leadingEdge_list[[x]],collapse = " "))) } 
    }
    colnames(leadingEdge)<-unique(DynamicPath %>% dplyr::select(DayType) %>% unlist())
    
    
    TypeUnique<-lapply(unique(tablePlot$type), function(x) paste(unique(tablePlot[tablePlot$type == x ,]$DayTypeSens), collapse = ".")) %>% unlist() 
    names(TypeUnique)<-unique(tablePlot$type)
    
    tablePlottest<-tablePlottest %>% mutate(type=factor(tablePlottest$type, levels = orderTypecell))
    col1<- tablePlot %>% arrange(type) %>% dplyr::select(color) %>% unlist() %>% as.character() %>% unique()

    listgenetoprint<-lapply(unique(tablePlottest$DayType), function(x) tablePlottest %>% filter(DayType == x) %>% arrange(padj,log2F)  %>% filter(padj <= 0.05) %>% dplyr::select(symbol)  %>% unique() %>% head(n=5) %>% unlist() )
    dt<-unique(tablePlottest$DayType)
    subset2<-data.frame()
    for(i in 1:length(dt) ){subset2<- rbind(subset2,tablePlottest %>% filter(DayType == dt[i]) %>% filter(symbol %in% listgenetoprint[[i]]))}
    col2<- lapply(sort(unique(subset2$type)) , function(x) colorsType[[x]]) %>% unlist() %>% as.character()

    mH_list<-matricePathHeatmap(path, tablePlottest)
    
    rmarkdown::render(paste0('/home/bioinfo/BulkAnalysis_plusNetwork/scripts/FigureSupp_gsea_Report.Rmd'),params = list("tablePlot"=DynamicPath,"tablePlottest"=tablePlottest, "subset"=subset,"TypeUnique"=TypeUnique,"col1"=col1,"col2"=col2, "minCountNormalized"=minCountNormalized,"maxCountNormalized"=maxCountNormalized, "minlogF"=minlogF,"maxlogF"=maxlogF,"title"=path, "subtitle"=NameDinamic,"numberGeneInPathway"=numberGeneInPathway,"tableLeadingEdge"=leadingEdge,"Heatmplot"=mH_list), output_file = paste0(dir,reportsuppPathway,"PathwayEnrichedIn",numberCond,"DayTypeCell_",path,"_meanCountNormalizedGene.html") )
    
    PlotGenelogFacrossDay <- ggplot(tablePlottest, aes(x=DayAge, y=meanCountNormalized,  color=log2F, group=GeneDayType)) +
      #scale_y_continuous(breaks = c(0,log2(10),log2(100),log2(1000),log2(10000)),labels=c("0",'10',"100","1000","10000")) + # fill=name allow to automatically dedicate a color for each group
      scale_y_log10()+
      geom_line(aes(col=log2F)) +
      geom_point(aes(shape=age,col=log2F))+
      scale_color_gradient2(midpoint=0,  low="blue", mid="grey88",high="red")+
      new_scale("color") +
      facet_wrap(~type,dir = "v", ncol=1, labeller = labeller(type=TypeUnique)) +
      geom_text_repel(data=subset2,
                      aes(label=symbol,color=type),  size=3,min.segment.length =0, segment.size = .8,point.padding= 0.8, force=5, force_pull = 0.1, max.overlaps=15,show.legend=F)  +
      scale_color_manual(values = col2, name="Top 10 genes padj<0.05")+
      ylab("meanCountNormalized")+ggtitle(unique(tablePlottest$Pathway))+
      theme_linedraw()+ theme(strip.text  = element_text(size = 12 , face = "bold"),strip.background = element_rect(fill="grey"),panel.grid.major.x = element_line(color = "gray60", size = 0.8),axis.title= element_text(face="bold"),axis.text = element_text(face="bold",size=14))
    
    g <- ggplot_gtable(ggplot_build(PlotGenelogFacrossDay))
    stripRowName <- which(grepl('strip', g$layout$name))
    k <- 1
    fills <- rev(col1)
    for (i in stripRowName) {
      j <- which(grepl('rect', g$grobs[[i]]$grobs[[1]]$childrenOrder))
      g$grobs[[i]]$grobs[[1]]$children[[j]]$gp$fill <- fills[k]
      k <- k+1
    }
    
    PlotGenelogFacrossDay2 <- ggplot(tablePlottest, aes(x=day, y=log2(meanCountNormalized+1), group=paste(symbol,type,age))) +
      scale_y_continuous(limits=c(minCountNormalized,maxCountNormalized),breaks = c(0,log2(10),log2(100),log2(1000),log2(10000)),labels=c("0",'10',"100","1000","10000")) + # fill=name allow to automatically dedicate a color for each group
      geom_line(aes(alpha=0.8) ,show.legend = F) +
      scale_color_gradient2(midpoint=0,  low="blue", mid="grey88",high="red", limits = c(minlogF,maxlogF),name="log2FoldChange")+
      geom_hline(yintercept = 0, color="black", linetype="dashed",lwd=0.5,slope=F)+
      geom_point(aes(col=log2F,size=log10(meanCountNormalized+1)))+
      scale_size(name="Counts normalised",limits=c(1,maxCountNormalized),breaks = c(log10(10),log10(100),log10(1000),log10(50000)),labels=c('10',"100","1000","50000"))+
      new_scale("color") +
      facet_grid(vars(type),vars(age),
                 labeller=labeller(type=TypeUnique))+
      geom_text_repel(data= subset2,
                      aes(label=symbol, color = type), size=3, segment.size = .1, force=5, force_pull = 5, max.overlaps=15,show.legend = F)  +
      scale_color_manual(values = col2, name="Top 10 genes padj<0.05")+
      ylab("meanCountNormalized")+ggtitle(unique(tablePlottest$Pathway))+
      theme_linedraw()+ theme(strip.text.x  = element_text(size = 12 , face = "bold"),strip.text.y  = element_text(size = 7 , face = "bold"),strip.background = element_rect(fill="grey"),panel.grid.major.x = element_line(color = "gray60", size = 0.8),axis.title= element_text(face="bold"),axis.text = element_text(face="bold",size=14))
    
    g2 <- ggplot_gtable(ggplot_build(PlotGenelogFacrossDay2))
    stripRowName <- which(grepl('strip', g2$layout$name))
    k <- 1
    fills <- c(rev(colorsAge),col1)
    for (i in stripRowName) {
      j <- which(grepl('rect', g2$grobs[[i]]$grobs[[1]]$childrenOrder))
      g2$grobs[[i]]$grobs[[1]]$children[[j]]$gp$fill <- fills[k]
      k <- k+1
    }
    
    save_plot(paste0(odir,plotpaways,"PathwayEnrichedIn_",numberCond,"_DayTypeCell_",path,"_meanCountNormalizedGene.png"),g,base_height = 10,base_width =14 )
    save_plot(paste0(odir,plotpaways,"PathwayEnrichedIn_",numberCond,"_DayTypeCell_",path,"_log2FGene.png"),g2,base_height = 10,base_width =14 )
    
    png(file=paste0(odir,plotpaways,"PathwayEnrichedIn_",numberCond,"_DayTypeCell_",path,"_heatmapScaleFixe.png"),units = "in", width=14, height= 9, res = 800, family = "Arial")
    draw (mH_list[[1]], column_title = path,
          heatmap_legend_side = "bottom", padding = unit(c(2, 2, 2, 6), "mm"))
    dev.off() 
    png(file=paste0(odir,plotpaways,"PathwayEnrichedIn_",numberCond,"_DayTypeCell_",path,"_heatmapScaleFree.png"),units = "in", width=14, height= 9, res = 800, family = "Arial")
    draw (mH_list[[2]], column_title = path,
          heatmap_legend_side = "bottom", padding = unit(c(2, 2, 2, 6), "mm"))
    dev.off() 
    }
  }
}

PathwaysMeanTMP8=PlotGeneCountNorm(Dynamic8table,msigdbr_list,8)
PathwaysMeanTMP7=PlotGeneCountNorm(Dynamic7table,msigdbr_list,7)
PathwaysMeanTMP6=PlotGeneCountNorm(Dynamic6table,msigdbr_list,6)
PathwaysMeanTMP5=PlotGeneCountNorm(Dynamic5table,msigdbr_list,5)
PathwaysMeanTMP4=PlotGeneCountNorm(Dynamic4table,msigdbr_list,4)
PathwaysMeanTMP3=PlotGeneCountNorm(Dynamic3table,msigdbr_list,3)
PathwaysMeanTMP2=PlotGeneCountNorm(Dynamic2table,msigdbr_list,2)
PathwaysMeanTMP1=PlotGeneCountNorm(Dynamic1table,msigdbr_list,1)


Dynamictable<-rbind(Dynamic8table,Dynamic7table,Dynamic6table,Dynamic5table,Dynamic4table,Dynamic3table,Dynamic2table,Dynamic1table)

leadingEdge_table<-data.frame()
pathway_list<-unique(Dynamictable$pathway)
gene_list<-lapply(1:length(pathway_list), function(x) msigdbr_list[pathway_list[x]] %>% unlist() %>% unique())
names(gene_list)<-pathway_list
DynamicTableSmall<- Dynamictable %>% dplyr::select(CondSigni,pathway,day,type,DayType,DayTypeSens,NES,color,leadingEdge)
for ( path in  pathway_list){
    tablePlot<-data_frame()
    numberGeneInPathway=length(unlist(gene_list[path]))
    DynamicPath<-DynamicTableSmall %>% filter(pathway==path) %>% mutate(type=factor(type,level=orderTypecell))
    NameDinamic<-unique(DynamicPath$CondSigni)
    tablePlotint<-cbind(symbol=unlist(gene_list[path]),do.call("rbind",replicate(numberGeneInPathway,DynamicPath,simplify = FALSE))%>%arrange(DayType))
    tablePlotint<-tablePlotint %>% mutate(GeneDayType=paste0(symbol,'_',day,"_",type))
    tablePlot<-rbind(tablePlot,tablePlotint)
    tablePlot<- tablePlot %>% mutate( padj = as.numeric(DEGtableTotalUPDOWN[match(tablePlot$GeneDayType,DEGtableTotalUPDOWN$GeneDayType),]$padj))
    tablePlot <- tablePlot %>% filter(padj < 0.05)
    
    geneDEGdt_list<-lapply(unique(DynamicPath %>% dplyr::select(DayType) %>% unlist()), function(x) tablePlot %>% filter(DayType == x) %>% dplyr::select(symbol) %>% unlist() %>% as.character())
    names(geneDEGdt_list)<-unique(DynamicPath %>%  dplyr::select(DayType) %>% unlist())
    leadingEdge_list <- lapply(unique(DynamicPath %>% dplyr::select(DayType) %>% unlist()), function(x) DynamicTableSmall %>% filter(pathway == path,DayType == x) %>% dplyr::select(leadingEdge) %>% unlist() %>% str_split(pattern = " " ) %>% unlist() %>% str_replace_all(pattern = "\n", replacement = "") %>% keep( ~ any(geneDEGdt_list[[x]] %in% .x)))  
    leadingEdge<-data.frame(x=c(length(leadingEdge_list[[1]]), paste(leadingEdge_list[[1]],collapse = " ")))
    
    leadingEdge_table<- rbind(leadingEdge_table, data.frame())
    if(length(leadingEdge_list) >1){
      for( x in 2:length(leadingEdge_list)){leadingEdge<-cbind(leadingEdge,c(length(leadingEdge_list[[x]]), paste(leadingEdge_list[[x]],collapse = " "))) } 
    }
    colnames(leadingEdge)<-unique(DynamicPath %>% dplyr::select(DayType) %>% unlist())
    
}
saveRDS(Dynamictable,paste0(odir,tableGSEA,TableGSEADynamics,'.rds'))
saveRDS(Dynamictable,paste0(odir,tableGSEA,TableGSEADynamics,'.csv'))

#####################################################
#### Looking pathways with hierarchical point of view
#####################################################

if (file.exists(tablePathwayHierarchieFull) == F ){

ConvertPathwaysName<-read.csv("data/ReactomePathways.txt",sep = "\t",header = F)
ConvertPathwaysName <- ConvertPathwaysName %>% filter(str_detect(V3,"Mus musculus"))

hierachicalPathway<-read.csv("data/ReactomePathwaysRelation.txt",sep = "\t",header = F)
hierachicalPathway<-hierachicalPathway%>% filter(str_detect(V1,"R-MMU")) %>% mutate("Node1"= str_replace_all(str_to_upper(ConvertPathwaysName[na.omit(match(hierachicalPathway$V1,ConvertPathwaysName$V1)),]$V2), "[\\ /\\-\\:]",'_'))  %>% mutate("Node1"= str_replace_all(Node1,"[,()]", ""))
hierachicalPathway<-hierachicalPathway%>% filter(str_detect(V2,"R-MMU")) %>% mutate("Node2"= str_replace_all(str_to_upper(ConvertPathwaysName[na.omit(match(hierachicalPathway$V2,ConvertPathwaysName$V1)),]$V2),"[\\ /\\-\\:]",'_')) %>% mutate("Node2"= str_replace_all(Node2,"[,()]", ""))

topHierarchie<-setdiff(hierachicalPathway$Node1,hierachicalPathway$Node2)
hierachicalPathway2<- rbind(cbind(hierachicalPathway$Node1,hierachicalPathway$Node2),cbind(rep("IN_REACTOME_HIERARCHIE",by=length(topHierarchie)),topHierarchie))

colnames(hierachicalPathway2)<-c("Node1","Node2")

GSEAtempo<-GSEAsigni%>% mutate( pathway2 = str_replace_all(pathway, "REACTOME_", "") )
listpathwaytotal<- c(unique(GSEAtempo$pathway2),"IN_REACTOME_HIERARCHIE")
listpathwayRestant<- listpathwaytotal

listpathcommun<-c()
for (path in listpathwaytotal){  if ( path %in% hierachicalPathway2[,2] ){ listpathcommun<-c(listpathcommun,path) }}
listPathwayNotInHierarchie<-setdiff(listpathwaytotal,listpathcommun)
listPathwayNotInHierarchie<-as.data.frame(listPathwayNotInHierarchie) %>% filter(!grepl("DISEASE|INFECTION|DEVELOPMENT|HIV|HCMV|CANCER|MEIOTIC|MEIOSIS|IN_REACTOME_HIERARCHIE",listPathwayNotInHierarchie))
NumberPathwayInHierarchie<-intersect(listpathwaytotal,listpathcommun)

if (length(listPathwayNotInHierarchie) > 0){
  #data_edged=as.data.frame(cbind(rep("root2", length(listPathwayNotInHierarchie$listPathwayNotInHierarchie)),listPathwayNotInHierarchie$listPathwayNotInHierarchie))
  hierachicalPathway2<- rbind(cbind(rep("NOT_IN_REACTOME_HIERARCHIE", length(listPathwayNotInHierarchie$listPathwayNotInHierarchie)),listPathwayNotInHierarchie$listPathwayNotInHierarchie),hierachicalPathway2)
  hierachicalPathway2<- rbind(c("root1","IN_REACTOME_HIERARCHIE"),hierachicalPathway2)
  hierachicalPathway2<- rbind(c("root1","NOT_IN_REACTOME_HIERARCHIE"),hierachicalPathway2)
}


explorer <- function(hierachicalPathway2,node2,path,listSommetparcouru,data_edged, level){
  level=level+1
  listSommetparcouru<-c(listSommetparcouru,node2)
  newpath=node2
  #size=GSEAtempo %>% filter(pathway==paste0("REACTOME_",node2)) %>% select(size) %>% unlist() %>% as.character() %>% as.numeric()
  size=ifelse(node2 == "root2", 1,GSEAtempo %>% filter(pathway==paste0("REACTOME_",node2)) %>% dplyr::select(size) %>% unlist() %>% as.character() %>% as.numeric() %>% sum(na.rm = T))
  data_edged<-rbind(data_edged,c(path,newpath,level,size))
  path=newpath
  for ( node2 in hierachicalPathway2[hierachicalPathway2[,1] == node2,2]){
    if (!(node2 %in% listSommetparcouru)){
      #print(c(node2,path))
      data_edged<-explorer(hierachicalPathway2,node2,path,listSommetparcouru,data_edged,level)
    }
  }
  return(data_edged)
}
data_edged_zoom=data.frame()
listSommetparcouru<-c()
data_edged_zoom<-explorer(hierachicalPathway2,"root1","root1",listSommetparcouru,data_edged_zoom,0)
colnames(data_edged_zoom)<-c("Node1","Node2","level","size")
write.table(data_edged_zoom,listHierarchieREACTOME)

data_edged_zoom2<-data_edged_zoom %>% filter(Node2!="root1") %>% dplyr::select(Node1,Node2)

data_tree <- FromDataFrameNetwork(data_edged_zoom2)
data_nested <- ToDataFrameTree(data_tree, 
                               level1 = function(x) x$path[2],
                               level2 = function(x) x$path[3],
                               level3 = function(x) x$path[4],
                               level4 = function(x) x$path[5],
                               level5 = function(x) x$path[6],
                               level6 = function(x) x$path[7],
                               level7 = function(x) x$path[8],
                               level8 = function(x) x$path[9],
                               level9 = function(x) x$path[10],
                               level10 = function(x) x$path[11],
                               level11 = function(x) x$path[12],
                               level12 = function(x) x$path[13],
                               level13 = function(x) x$path[14])[-1,-1]
# Now we can plot it as seen before!
data_nested<-data_nested %>% arrange(desc(level13),desc(level12),desc(level11),desc(level10),desc(level9),desc(level8),desc(level7),desc(level6),desc(level5),desc(level4),desc(level3),desc(level2),desc(level1))

listpathway<-listpathwaytotal

data_nested_print<-data.frame()
data_nested_full<-data.frame()
for(l in rownames(data_nested)){
  newrow=data_nested[l,]
  newrowFull=data_nested[l,]
  pathTF=F
  for (c in 13:1){
    if (is.na(data_nested[l,c]) == F ){
      if ( data_nested[l,c] %in% listpathway){
        pathTF=T
        newrow[[c]]<- paste0('<b>',str_replace_all(data_nested[l,c],'_',' '), "</b>")
        }
      else{
        if (pathTF == F){
            newrow[[c]]<- ""
            newrowFull[[c]]<-""
        }
        else{
          newrow[[c]]<- tolower(str_replace_all(data_nested[l,c],'_',' '))
        }
      }
    }
  }
  if (pathTF == T){
    data_nested_print<-rbind(data_nested_print,newrow)
    data_nested_full<-rbind(data_nested_full,newrowFull)
  }
}
data_nested_print<-data_nested_print %>% arrange(desc(level13),desc(level12),desc(level11),desc(level10),desc(level9),desc(level8),desc(level7),desc(level6),desc(level5),desc(level4),desc(level3),desc(level2),desc(level1))
data_nested_full<-data_nested_full %>% arrange(desc(level13),desc(level12),desc(level11),desc(level10),desc(level9),desc(level8),desc(level7),desc(level6),desc(level5),desc(level4),desc(level3),desc(level2),desc(level1))


write.table(data_nested_full, tablePathwayHierarchieFull )
write.table(data_nested_print, tablePathwayHierarchieTotalPathwayEnrichie )
}else{
  data_nested_full=read.table(tablePathwayHierarchieFull )
  data_nested_print=read.table(tablePathwayHierarchieTotalPathwayEnrichie )
  data_edged_zoom=read.table(listHierarchieREACTOME)
}

data_nested_print[is.na(data_nested_print)] <- ""

pathwayCheck<-c()
rowsupp<-c()
for(l in rownames(data_nested_print)){
  pathTF=F
  for (c in 13:1){
    if (data_nested_print[l,c]!= ""){
      if (pathTF == F){
        if (data_nested_print[l,c] %in% pathwayCheck){
          rowsupp<-c(rowsupp,l)
          break
        }
      }
      pathTF=T
      pathwayCheck<-c(pathwayCheck,data_nested_print[l,c])
    }
  }
}
data_nested_3<-data_nested_print[-as.numeric(rowsupp),]
data_nested_3<-data_nested_3 %>% arrange(level1,level2,level3,level4,level5,level6,level7,level8,level9,level10,level11,level12,level13)
write.table(data_nested_3, file=paste0( odir,HierarchieData,"ALL_table_REACTOME_Hierachie_print.txt"))

data_nested_full[is.na(data_nested_full)] <- ""


for (t in Typecellv){
  for (d in days){
  #prefix_name<-paste0(t)
  #prefix_name<-paste0(d)
  prefix_name<-paste0(d,"_",t)
  #prefix_name<-paste0("ALL")
  
  
  #GSEAtempo<-GSEAsigni %>% filter(type==t) %>% mutate( pathway2 = str_replace_all(pathway, "REACTOME_", "") )
  #GSEAtempo<-GSEAsigni %>% filter(day==d) %>% mutate( pathway2 = str_replace_all(pathway, "REACTOME_", "") )  
  GSEAtempo<-GSEAsigni%>% filter(type==t,day==d) %>% mutate( pathway2 = str_replace_all(pathway, "REACTOME_", "") )
  #GSEAtempo<-GSEAsigni %>% mutate( pathway2 = str_replace_all(pathway, "REACTOME_", "") )
  
  if (dim(GSEAtempo)[1] >0){
    listpathway<- c(unique(GSEAtempo$pathway2),"NOT_IN_REACTOME_HIERARCHIE","IN_REACTOME_HIERARCHIE")
    data_nested_2<-data.frame()
    sizeVector<-c()
    for( l in rownames(data_nested_full)){
      newrow=rep("", 13)
      sizetempo=0
      pathFound=F
      l=as.numeric(l)
      for (c in 13:1){
        if (data_nested_full[l,c]!= ""){
          if ( data_nested_full[l,c] %in% listpathway){
            if (pathFound == F){
              sizetempo=GSEAtempo %>% filter(pathway2 == data_nested_full[l,c] ) %>% dplyr::select( size ) %>% unlist() %>% as.character() %>% as.numeric()
              sizetempo=sum(sizetempo,na.rm = T )
              }
            pathFound=T 
            newrow[[c]]<- str_replace_all(data_nested_full[l,c],'_',' ')
          }else{
            if (pathFound == T){
              newrow[[c]]<- tolower(str_replace_all(data_nested_full[l,c],'_',' '))
            }
          }
        }
        }
      if (pathFound == T){
        newrowExist<-t(data.frame(newrow))
        colnames(newrowExist)<-colnames(data_nested_2)
        if (nrow(merge(newrowExist,data_nested_2))==0){
          data_nested_2<-rbind(data_nested_2,newrow)
          sizeVector<-c(sizeVector,unique(sizetempo))
        }
      }
    }

    colnames(data_nested_2)<-c("level1","level2","level3","level4","level5","level6","level7","level8","level9","level10","level11","level12","level13")
    data_nested_2$value=sizeVector
    
    if (dir.exists(paste0(odir,HierarchieData,prefix_name,"/")) == F) {
      dir.create(paste0(odir,HierarchieData,prefix_name,"/"))
    }
    write.table(data_nested_2, file=paste0( odir,HierarchieData,prefix_name,"/",prefix_name, "_table_REACTOME_Hierachie.txt"))

    data_nested_2$pathString <- paste("roots", data_nested_2$level1, data_nested_2$level2, data_nested_2$level3, data_nested_2$level4,data_nested_2$level5,data_nested_2$level6,data_nested_2$level7, data_nested_2$level8, data_nested_2$level9, data_nested_2$level10,data_nested_2$level11,data_nested_2$level12, sep = "/")
    data_Node <-as.Node(data_nested_2)
    if (t ==""){
      p <- circlepackeR(data_Node, size = "value",color_min = paste0("hsl(56,80%,80%)"), color_max = paste0("hsl(341,30%,30%)"), width = 1800, height = 1000)
    }else{
    p <- circlepackeR(data_Node, size = "value",color_min = paste0("hsl(",hsl_color[[t]],",80%,80%)"), color_max = paste0("hsl(",hsl_color[[t]],",30%,30%)"), width = 1500, height = 1000)
    }
    p
    saveWidget(p, file=paste0( odir,plotpathwaysHierarchie,"/",prefix_name, "_REACTOME_Hierachie.html"))
  }
}
}
for (t in Typecellv){
  for(d in daysv[[t]]){
    prefix_name<-paste0(d,"_",t)
    GSEAtempo<-GSEAsigni%>% filter(day==d,type==t) %>% mutate( pathway2 = str_replace_all(pathway, "REACTOME_", "") )
    if (dim(GSEAtempo)[1] >0){
      rmarkdown::render(paste0('/home/bioinfo/BulkAnalysis_plusNetwork/scripts/FigureSupp_gsea_Hierarchie.Rmd'),params = list("title"=paste0("Reactome_hierarchie_in_",prefix_name),"prefix_name"=prefix_name,"GSEAsigni"=GSEAtempo), output_file = paste0(dir,reportsuppHierarchie,"/",prefix_name,"_Hierarchical_Figsupp") )
    }
  }
}

for (t in Typecellv){
    prefix_name<-paste0(t)
    GSEAtempo<-GSEAsigni%>% filter(type==t) %>% mutate( pathway2 = str_replace_all(pathway, "REACTOME_", "") )
    if (dim(GSEAtempo)[1] >0){
      rmarkdown::render(paste0('/home/bioinfo/BulkAnalysis_plusNetwork/scripts/FigureSupp_gsea_Hierarchie.Rmd'),params = list("title"=paste0("Reactome_hierarchie_in_",prefix_name),"prefix_name"=prefix_name,"GSEAsigni"=GSEAtempo), output_file = paste0(dir,reportsuppHierarchie,"/",prefix_name,"_Hierarchical_Figsupp") )
  }
}

for(d in days){
    prefix_name<-paste0(d)
    GSEAtempo<-GSEAsigni%>% filter(day==d) %>% mutate( pathway2 = str_replace_all(pathway, "REACTOME_", "") )
    if (dim(GSEAtempo)[1] >0){
      rmarkdown::render(paste0('/home/bioinfo/BulkAnalysis_plusNetwork/scripts/FigureSupp_gsea_Hierarchie.Rmd'),params = list("title"=paste0("Reactome_hierarchie_in_",prefix_name),"prefix_name"=prefix_name,"GSEAsigni"=GSEAtempo), output_file = paste0(dir,reportsuppHierarchie,"/",prefix_name,"_Hierarchical_Figsupp") )
    }
  }
  



rmarkdown::render("~/BulkAnalysis_plusNetwork/scripts/fgseaReport.Rmd", output_dir = "~/BulkAnalysis_plusNetwork/reports/")

