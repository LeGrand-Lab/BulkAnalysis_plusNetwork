# Organize GSEA results into tables
#   check fgsea_dayandcelltype.R
# output dir: exam_INTER_conditions/static/GSEA/ 
# and pdf plots
# --
# johaGL
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

##################################
#### Path files to load or created
##################################

#Directories
setwd("~/BulkAnalysis_plusNetwork2/")
analyse<-"analysis/AnalysisGSEA/"
odir <- "/home/bioinfo/BulkAnalysis_plusNetwork2/exam_INTER_conditions/static/"
gseardsFull = "GSEA/rds/"
gseacsvFull= "GSEA/csv/"
HierarchieData= "GSEA/HierarchieData/"

#Loads
gseaoutfull2padj = "fgsea_padj_full.rds" 
gseaout_filtered_full2padj = "fgsea_padj_filtered.rds" 

#Created
GSEAconcatPath<-paste0(odir,gseardsFull,"GSEAconcat.rds")
GSEAconcatsigniPath<-paste0(odir,gseardsFull,"GSEAconcatSigni.rds")
CountNormalizedYoungPath<-paste0(odir,"csv/YoungCountNormalized")
CountNormalizedOldPath<-paste0(odir,"csv/OldCountNormalized")
listHierarchieREACTOME<-paste0(odir,HierarchieData,"ListReactomeHiearachieUseFunction.txt")
tablePathwayHierarchieFull<-paste0(odir,HierarchieData,"TableReactomeHiearachieUseFunction.txt")
tablePathwayHierarchieTotalPathwayEnrichie<-paste0(odir,HierarchieData,"TableReactomeHiearachieTotalPathwayEnrichie.txt")

if (dir.exists(paste0(odir,analyse,"PathwayEnrichedIn",numberCond,"DayTypeCell/")) == F) {
  dir.create(paste0(odir,analyse,"PathwayEnrichedIn",numberCond,"DayTypeCell/"))
}

##################################
#### Load files
##################################

tablesneeded <- F
mat4heatmapneeded <- F
fullrds <- readRDS(paste0(odir,gseardsFull, gseaoutfull2padj))
pathsFiltered = readRDS( paste0(odir, gseardsFull,  gseaout_filtered_full2padj) )


thegmt <- msigdbr(species = "Mus musculus", 
                  category = 'C2', 
                  subcategory=c('CP:REACTOME'))
msigdbr_list = split(x = thegmt$gene_symbol, f = thegmt$gs_name)

fullDEsta = readRDS(paste0(odir, "rds/shot_rds_full.rds"))

genes_df <- read.table("data/genesinfo.csv", sep="\t", header=T)

daysv<-list()
daysv<-lapply( unique(fullDEsta$type), function(t) unique(fullDEsta[fullDEsta$type == t,]$day) )
names(daysv)<-unique(fullDEsta$type)
Typecellv<-unique(fullDEsta$type)
days<-unique(fullDEsta$day)

####################
####Function colors
####################
darker.col = function(color, how.much = 30){
  colorRampPalette(c("white",color,"black"))(9)[how.much]
}
tablecolor<-data.frame("ECs"=darker.col("#0072B2", c(4:7,5)),"FAPs"=darker.col("#F0E442", c(3:6,7)),"M1"=darker.col("#D55E00", c(3:6,6)),"M2"=darker.col("#CC79A7", c(3:6,6)),"Neutro"=darker.col("#009E73", c(3:6,6)),"sCs"=darker.col("#56B4E9",c(2:5,6)))
rownames(tablecolor)<-c("D0","D2","D4","D7","Police")
hsl_color<-data.frame("ECs"=270,"FAPs"=56,"M1"=26,"M2"=327,"Neutro"=164,"sCs"=202)


###################################################
####Extract GSEA individual table and merge in xlsx
####################################################

if (tablesneeded){
  system(paste0("cd ",odir,"GSEA/; ", 
               "if [ ! -d csv/ ]; then mkdir csv; else echo 'csv/ exists, no action';fi"))
  days = names(fullrds)
  OUT <- createWorkbook()
  for (d in names(fullrds)){
    print(paste("\n^^^^^", d, "^^^^^"))
    for (t in names(fullrds[[d]])){
      tmpdf <- fullrds[[d]][[t]]
      tmpdf <- tmpdf %>% 
        mutate(gene_symbols = sapply(leadingEdge, function(x) paste(x, collapse=", ")),
               .before = leadingEdge) %>%
        select(-leadingEdge)
      tmpdf <- tmpdf %>% arrange(padj,desc(abs(NES)),desc(length(gene_symbols)))
      
      write.table(tmpdf, paste0(odir,gseacsvFull, t,"_",d,"_pathways.csv"), sep="\t", 
                  col.names = T, row.names=F              )
      addWorksheet(OUT, paste0(t,'_',d))
      writeData(OUT,sheet= paste0(t,'_',d),x=tmpdf)
    }
  }
  saveWorkbook(OUT,paste0(odir,gseacsvFull,"Total_GSE_pathways.xlsx"),overwrite = T)
  rm(fullrds)
}

####################################################
# Table GSEA concat to explore combinaison pathway
####################################################

if (file.exists(GSEAconcatsigniPath) == F ){

XX_full = readRDS(paste0(odir,gseardsFull, gseaoutfull2padj))
GSEAconcat<-tibble() # table concat to explore combinaison pathway signi
for ( t in Typecellv ){
  for (d in daysv[[t]]){
    GSEAconcat<-rbind(GSEAconcat, XX_full[[d]][[t]] %>% mutate(day=d) %>% mutate(type=t))
  }
}

GSEAsigni<- GSEAconcat %>% filter(padj<=0.05) %>% arrange(day,type) %>% mutate(sens=ifelse(NES < 0,"down", "up"))

saveRDS(GSEAconcat,GSEAconcatPath)
saveRDS(GSEAsigni,GSEAconcatsigniPath)

}else{
  GSEAconcat=readRDS(GSEAconcatPath)
  GSEAsigni=readRDS(GSEAconcatsigniPath)
}


################################################# 
# Extract stats significant pathways enrichied
#################################################
NumberPathwaysSigni = array(NA, dim=c(length(days),length(Typecellv)))
rownames(NumberPathwaysSigni) =  days
colnames(NumberPathwaysSigni) = Typecellv
NumberPathwaysUniqueOnOneDay<-c()
NumberPathwaysUniqueOnOneTypeCell<-c()
for (d in days){
  NumberPathwaysUniqueOnOneDay<-c(NumberPathwaysUniqueOnOneDay,GSEAsigni %>% filter(day==d) %>% select(pathway) %>% unique() %>% unlist() %>% length())
  for (t in Typecellv){
    NumberPathwaysSigni[d, t] <- dim(GSEAsigni %>% filter(day==d & type==t))[1]
    if (d == days[1]){
      NumberPathwaysUniqueOnOneTypeCell<-c(NumberPathwaysUniqueOnOneTypeCell, GSEAsigni %>% filter(type==t) %>% select(pathway) %>% unique() %>% unlist() %>% length())
    }
  }
}

names(NumberPathwaysUniqueOnOneDay)<-days
names(NumberPathwaysUniqueOnOneTypeCell)<-Typecellv

#### Explore Dynamique pathways Enrichied padj<0.05

## table of number of pathways found in 1, 2 , ... conditions
RepartitionPathwaySummary<-table(table(GSEAsigni$pathway))
## dataframe : Var1 = pathways Freq = Number condition which the pathway is found
RepartitionPathway<-as.data.frame(table(GSEAsigni$pathway)) 

## NewTable with NbCondSigni, CondSigni: Day1_Type1_Sens1.Day_Type2_Sens2 ..., DayTypeSens,NbPathWithDynamique,color
DynamicPathwayAcrossCond <- function(RepartitionPathway,GSEAsigni,n){
  PathwaySigniCond<-RepartitionPathway %>% filter(Freq==n)
  PathwaySigniCondtable<- lapply(PathwaySigniCond$Var1 , function(x) dplyr::filter(GSEAsigni, pathway == x))
  DynamicCond <-lapply(PathwaySigniCondtable, function(x) paste0(x$day,'_',x$type,'_',x$sens) %>% sort())
  names(DynamicCond)<-lapply(PathwaySigniCondtable, function(x) unique(x$pathway) )
  DynamicCondconcat<-c()
  for(Dym in DynamicCond){DynamicCondconcat=c(DynamicCondconcat,str_c(Dym,collapse = " "))}
  names(DynamicCondconcat)<-lapply(PathwaySigniCondtable, function(x) unique(x$pathway))
  Dynamiquetable<-ldply(lapply(PathwaySigniCondtable, function(x) as.data.frame(cbind(rep(n,n),x,paste0(x$day,'_',x$type,'_',x$sens),rep(DynamicCondconcat[[unique(x$pathway)]],n)))),rbind)
  colnames(Dynamiquetable)<-c("NbCondSigni",colnames(PathwaySigniCondtable[[1]]),"DayTypeSens","CondSigni")
  
  Dynamique<-lapply(unique(Dynamiquetable$pathway), function(x) unique(Dynamiquetable %>% filter(pathway==x) %>% select(CondSigni) %>% unlist())) %>% unlist()
  names(Dynamique)<-unique(Dynamiquetable$pathway)
  tableDynamique<-table(Dynamique)
  Dynamiquetable <- Dynamiquetable %>% mutate(NbPathWithDynamique=as.numeric(as.character(tableDynamique[match(Dynamiquetable$CondSigni,names(tableDynamique))]))) %>% mutate(color=as.character(tablecolor["D4",match(Dynamiquetable$type,colnames(tablecolor))])) %>% arrange(desc(NbPathWithDynamique))
  Dynamiquetable <- Dynamiquetable %>% mutate( DayType = paste0(day,'_',type))
  Dynamiquetable <-  Dynamiquetable %>% mutate( link=paste0("<a href='","../exam_INTER_conditions/static/",analyse,"PathwayEnrichedIn",n,"DayTypeCell/",pathway,"_meanCountNormalizedGene.html","'>FigureSupp</a>"))
  Dynamiquetable <-  Dynamiquetable %>% mutate( NamePathway= str_replace(str_replace_all(pathway,"_"," "), "REACTOME",""))
  
    return(Dynamiquetable)
}


#Plot of Pathway in function NES
GetPlot<-function(Dynamic6table,GSEAsigni, nbCond , vertical){
  col<- Dynamic6table %>% arrange(type) %>% select(color) %>% unlist() %>% as.character() %>% unique()
  Dynamic6table <-Dynamic6table %>% mutate( CondSigni = as.factor(CondSigni)) %>% mutate(CondSigni = fct_reorder(CondSigni, as.numeric(NbPathWithDynamique), .desc = T) )
  
  typebytime<-ggplot(Dynamic6table, aes(x=day,y=as.numeric(NES),fill=type,color=type,group=paste0(pathway,'_',type))) +
    geom_point(aes(size=as.numeric(NbPathWithDynamique),alpha=0.7),show.legend = c(size=T,alpha=F))+
    scale_size_binned(name="Number path",limits=c(1,MaxNbPaathWithSameDyn))+
    ylab("NES")+
    geom_line()+
    geom_hline(yintercept = 0, color="black", linetype="dashed",lwd=0.5,slope=F)+
    scale_color_manual(values=col)+
    ggtitle( paste0("Pathways with ",nbCond," conditions enriched")) +
    facet_wrap(~CondSigni+NbPathWithDynamique ,dir = "v", ncol=vertical,labeller = labeller(NbPathWithDynamique = label_facet(Dynamic6table$NbPathWithDynamique,"Number pathways: ")))+
    theme_linedraw()+ theme(strip.text  = element_text(size = 11 , face = "bold"),strip.background = element_rect(fill="grey"),panel.grid.major.x = element_line(color = "gray60", size = 0.8),axis.title= element_text(face="bold"),axis.text = element_text(face="bold",size=14))
  print(typebytime)
  save_plot(paste0(odir, analyse, 'Pathway_enriched_in__',nbCond ,'conditions_Plot2.svg'),typebytime, base_height = 10,
            base_width = 14)
}

Dynamic8table=DynamicPathwayAcrossCond(RepartitionPathway,GSEAsigni,8)
Dynamic7table=DynamicPathwayAcrossCond(RepartitionPathway,GSEAsigni,7)
Dynamic6table=DynamicPathwayAcrossCond(RepartitionPathway,GSEAsigni,6)
Dynamic5table=DynamicPathwayAcrossCond(RepartitionPathway,GSEAsigni,5)
Dynamic4table=DynamicPathwayAcrossCond(RepartitionPathway,GSEAsigni,4)
Dynamic3table=DynamicPathwayAcrossCond(RepartitionPathway,GSEAsigni,3)
Dynamic2table=DynamicPathwayAcrossCond(RepartitionPathway,GSEAsigni,2)
Dynamic1table=DynamicPathwayAcrossCond(RepartitionPathway,GSEAsigni,1)

MaxNbPaathWithSameDyn<-max(Dynamic8table$NbPathWithDynamique,Dynamic7table$NbPathWithDynamique,Dynamic6table$NbPathWithDynamique,Dynamic5table$NbPathWithDynamique,Dynamic4table$NbPathWithDynamique,Dynamic3table$NbPathWithDynamique,Dynamic2table$NbPathWithDynamique,Dynamic1table$NbPathWithDynamique )

GetPlot(Dynamic8table,GSEAsigni, 8 , 1)
GetPlot(Dynamic7table,GSEAsigni, 7 , 1)
GetPlot(Dynamic6table,GSEAsigni, 6 , 1)
GetPlot(Dynamic5table,GSEAsigni, 5 , 2)
GetPlot(Dynamic4table,GSEAconcat, 4 , 2)
GetPlot(Dynamic3table,GSEAsigni, 3 , NA)
GetPlot(Dynamic2table,GSEAsigni, 2 , NA)
GetPlot(Dynamic1table,GSEAsigni, 1 , NA)

##############
# Explore data Expression associated to dynamic tendance found in enriched pathway
##############
if ( file.exists(CountNormalizedOldPath) == F){
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
write.table(CountNormalizedYoung,CountNormalizedYoungPath)
write.table(CountNormalizedOld,CountNormalizedOldPath)
print("created")
}else{
  print("load")
  CountNormalizedYoung=read.table(CountNormalizedYoungPath)
  CountNormalizedOld=read.table(CountNormalizedOldPath)
}

minCountNormalized<-log2(round(min(as.numeric(CountNormalizedYoung$countNormalized,CountNormalizedOld$countNormalized, na.rm=T )),digits = 0)+1)
maxCountNormalized<-log2(round(max(as.numeric(CountNormalizedYoung$countNormalized,CountNormalizedOld$countNormalized, na.rm=T )),digits = 0)+1)
  
minlogF <-round(min(fullDEsta$log2FoldChange, na.rm=T ),digits = 0)
maxlogF <- round(max(fullDEsta$log2FoldChange, na.rm=T ),digits = 0)

fullDEsta<- fullDEsta %>% unite( 'GeneDayType' , c(symbol,day,type) ,sep = '_',remove = F)


matricePathHeatmap<-function(path, tablePlottest){
  genes<-data.frame(symbol=NA)
  ntop=15
  while (length(genes$symbol) <= 40 ){
    if( length(genes$symbol) >= length(unique(tablePlottest$symbol)) )
      break
  genes<-lapply(unique(tablePlottest$DayType), function(x) tablePlottest %>% filter(DayType == x) %>% arrange(padj,log2F) %>% head(n=ntop) ) %>% purrr::reduce(full_join) %>% filter(padj <= 0.05) %>% select(symbol) %>% unique()
  ntop=ntop+1
  }
  matrix<-data.frame(symbol=genes$symbol )
  colNAmesMatrix<-c("symbol")
  typev<-c()
  for ( t in Typecellv){
    for ( d in daysv[[t]]){
      typev<-c(typev,t)
      colNAmesMatrix<-c(colNAmesMatrix,paste0(t,'_',d))
      if (length(fullDEsta %>% filter(day==d,type==t) %>% subset(symbol %in% genes$symbol) %>% select(log2FoldChange) %>% unlist() %>% as.numeric()) != length(genes$symbol)){
        tempo=fullDEsta %>% filter(day==d,type==t) %>% subset(symbol %in% genes$symbol)
        fullDEsta.tempo <- data.frame(tempo ,genes_df[match(tempo$id, genes_df$Geneid),]$symbol)
        extractGeneDuplicate <- fullDEsta.tempo[duplicated(tempo$symbol),]$symbol
        genes2 <-genes[! genes$symbol == extractGeneDuplicate,]
        
        genes<-data_frame("symbol"=genes2) %>% as.data.frame()
        
        matrix<-data.frame(symbol=genes$symbol )
        fullDEsta %>% filter(day==d,type==t) %>% subset(symbol %in%  genes$symbol) %>% select(symbol) %>% unlist()
        
        matrix <- matrix %>% mutate("random" = fullDEsta %>% filter(day==d,type==t) %>% subset(symbol %in%  genes$symbol) %>% select(log2FoldChange) %>% unlist() %>% as.numeric())
      }
      else{
      matrix <- matrix %>% mutate("random" = fullDEsta %>% filter(day==d,type==t) %>% subset(symbol %in% genes$symbol) %>% select(log2FoldChange) %>% unlist() %>% as.numeric() )
      }
      colnames(matrix)<-colNAmesMatrix
    }
  }
  rownames(matrix)<-genes$symbol
  matrix<-matrix[,-1]
  colgroup<-as.data.frame(typev)
  rownames(colgroup)<-colNAmesMatrix[2:length(colNAmesMatrix)]
  ha = HeatmapAnnotation(foo = anno_block(gp = gpar(fill = c("darkblue", "orange",
                                                             "firebrick","violet",
                                                             "darkgreen","royalblue") )))
  heatmap<-ComplexHeatmap::Heatmap(as.matrix(matrix),cluster_columns = FALSE,column_split = colgroup,top_annotation = ha,column_names_rot = 45,name=path,heatmap_legend_param = list(title="Log2FoldChange",
                                                                                                                                                                            direction = "horizontal"),use_raster = TRUE,col=colorRamp2(c(minlogF,0,maxlogF), c("blue", "white", "red")))
  return(heatmap)
  }

PlotGeneCountNorm <- function(Dynamic7table,msigdbr_list,numberCond){
  pathway_list<-unique(Dynamic7table$pathway)
  gene_list<-lapply(1:length(pathway_list), function(x) msigdbr_list[pathway_list[x]] %>% unlist() %>% unique())
  names(gene_list)<-pathway_list
  DynamicTableSmall<- Dynamic7table %>% select(CondSigni,pathway,day,type,DayType,DayTypeSens,NES,color)
  for ( path in  pathway_list){
    tablePlot<-data_frame()
    numberGeneInPathway=length(unlist(gene_list[path]))
    DynamicPath<-DynamicTableSmall %>% filter(pathway==path)
    NameDinamic<-unique(DynamicPath$CondSigni)
    tablePlotint<-cbind(symbol=unlist(gene_list[path]),do.call("rbind",replicate(numberGeneInPathway,DynamicPath,simplify = FALSE))%>%arrange(DayType))
    tablePlotint<-tablePlotint %>% mutate(GeneDayType=paste0(symbol,'_',day,"_",type))
    tablePlot<-rbind(tablePlot,tablePlotint)
    #fullDEsta<- fullDEsta %>% unite( 'GeneDayType' , c(symbol,day,type) ,sep = '_',remove = F)
    tablePlot<- tablePlot %>% mutate( padj = as.numeric(fullDEsta[match(tablePlot$GeneDayType,fullDEsta$GeneDayType),]$padj))
    tablePlot<- tablePlot %>% mutate( log2F = as.numeric(fullDEsta[match(tablePlot$GeneDayType,fullDEsta$GeneDayType),]$log2FoldChange))
    tablePlotY<- tablePlot %>% mutate( meanCountNormalized = as.numeric(CountNormalizedYoung[match(tablePlot$GeneDayType,CountNormalizedYoung$GeneDayType),]$countNormalized)) %>% mutate( age = rep("Young",length(tablePlot$day)))
    tablePlotO<- tablePlot %>% mutate( meanCountNormalized = as.numeric(CountNormalizedOld[match(tablePlot$GeneDayType,CountNormalizedOld$GeneDayType),]$countNormalized)) %>% mutate( age = rep("Old",length(tablePlot$day)))
    tablePlot <- rbind(tablePlotY,tablePlotO) %>% filter(padj < 0.05)
    tablePlottest <- tablePlot %>% unite("DayAge",c(day, age),sep = ' ',remove = F) 
    tablePlottest <- tablePlottest %>% arrange(day,desc(age)) 
    tablePlottest <- tablePlottest %>% mutate (DayAge = factor( tablePlottest$DayAge , levels =unique(tablePlottest$DayAge) ))
    
    leadingEdge <- lapply(unique(tablePlottest %>% arrange(type) %>% select(DayType) %>% unlist()), function(x) Dynamic7table %>% filter(pathway == path,DayType == x) %>% select(leadingEdge)) %>% as.data.frame()  
    names(leadingEdge)<-unique(tablePlottest %>% arrange(type) %>% select(DayType) %>% unlist())
    leadingEdge <- rbind(lapply(names(leadingEdge), function(x) length(leadingEdge[,x][[1]] )) %>% unlist(),leadingEdge)
    
    
    TypeUnique<-lapply(unique(tablePlot$type), function(x) paste(unique(tablePlot[tablePlot$type == x ,]$DayTypeSens), collapse = ".")) %>% unlist() 
    names(TypeUnique)<-unique(tablePlot$type)
    
    col1<- tablePlot %>% arrange(type) %>% select(color) %>% unlist() %>% as.character() %>% unique()
    
    subset<-lapply(unique(tablePlottest$DayType), function(x) tablePlottest %>% filter(DayType == x) %>% arrange(padj,log2F) %>% head(n=10L) ) %>% purrr::reduce(full_join) %>% filter(padj <= 0.05)
    col2<- lapply(sort(unique(subset$type)) , function(x) tablecolor["Police",x]) %>% unlist() %>% as.character()
    
    print(col1)
    print(col2)
    mH<-matricePathHeatmap(path, tablePlottest)
    
    rmarkdown::render(paste0('/home/bioinfo/BulkAnalysis_plusNetwork2/scripts/FigureSupp_gsea_Report.Rmd'),params = list("tablePlot"=DynamicPath,"tablePlottest"=tablePlottest, "subset"=subset,"TypeUnique"=TypeUnique,"col1"=col1,"col2"=col2, "minCountNormalized"=minCountNormalized,"maxCountNormalized"=maxCountNormalized, "minlogF"=minlogF,"maxlogF"=maxlogF,"title"=path, "subtitle"=NameDinamic,"numberGeneInPathway"=numberGeneInPathway,"tableLeadingEdge"=leadingEdge,"Heatmplot"=mH), output_file = paste0(odir,analyse,"PathwayEnrichedIn",numberCond,"DayTypeCell/",path,"_meanCountNormalizedGene.html") )
    
    PlotGenelogFacrossDay <- ggplot(tablePlottest, aes(x=DayAge, y=log2(meanCountNormalized+1),  color=log2F, group=GeneDayType)) +
      scale_y_continuous(limits=c(minCountNormalized,maxCountNormalized),breaks = c(0,log2(10),log2(100),log2(1000),log2(10000)),labels=c("0",'10',"100","1000","10000")) + # fill=name allow to automatically dedicate a color for each group
      geom_line(aes(col=log2F)) +
      geom_point(aes(shape=age,col=log2F))+
      scale_color_gradient2(midpoint=0,  low="blue", mid="grey88",high="red", limits = c(minlogF,maxlogF))+
      new_scale("color") +
      geom_text_repel(data= subset,
                      aes(label=symbol, color = type), size=3, segment.size = .1, force=2, force_pull = 2, max.overlaps=15)  +
      facet_wrap(~type,dir = "v", ncol=1, labeller = labeller(type=TypeUnique)) +
      scale_color_manual(values = col2, name="Top 10 genes padj<0.05")+
      ylab("meanCountNormalized")+ggtitle(unique(tablePlottest$Pathway))+
      theme_linedraw()+ theme(strip.text  = element_text(size = 12 , face = "bold"),strip.background = element_rect(fill="grey"),panel.grid.major.x = element_line(color = "gray60", size = 0.8),axis.title= element_text(face="bold"),axis.text = element_text(face="bold",size=14))
    print(PlotGenelogFacrossDay)
    
    PlotGenelogFacrossDay2 <- ggplot(tablePlottest, aes(x=day, y=log2F, group=paste(symbol,type,age))) +
      scale_y_continuous(limits=c(minlogF,maxlogF),) + # fill=name allow to automatically dedicate a color for each group
      geom_line() +
      scale_color_gradient2(midpoint=0,  low="blue", mid="grey88",high="red", limits = c(minlogF,maxlogF),name="log2FoldChange")+
      geom_hline(yintercept = 0, color="black", linetype="dashed",lwd=0.5,slope=F)+
      geom_point(aes(col=log2F,size=log10(meanCountNormalized+1)))+
      scale_size(name="Counts normalised",limits=c(1,maxCountNormalized),breaks = c(log10(10),log10(100),log10(1000),log10(50000)),labels=c('10',"100","1000","50000"))+
      new_scale("color") +
      geom_text_repel(data= subset,
                      aes(label=symbol, color = type), size=3, segment.size = .1, force=5, force_pull = 5, max.overlaps=15)  +
      facet_grid(vars(type),vars(age),
                 labeller=labeller(type=TypeUnique))+
      scale_color_manual(values = col2, name="Top 10 genes padj<0.05")+
      ylab("log2FoldChange")+ggtitle(unique(tablePlottest$Pathway))+
      theme_linedraw()+ theme(strip.text  = element_text(size = 12 , face = "bold"),strip.background = element_rect(fill="grey"),panel.grid.major.x = element_line(color = "gray60", size = 0.8),axis.title= element_text(face="bold"),axis.text = element_text(face="bold",size=14))
    print(PlotGenelogFacrossDay2)
    
    save_plot(paste0(odir,analyse,"PathwayEnrichedIn",numberCond,"DayTypeCell/",path,"_meanCountNormalizedGene.png"),PlotGenelogFacrossDay,base_height = 10,base_width =14 )
    save_plot(paste0(odir,analyse,"PathwayEnrichedIn",numberCond,"DayTypeCell/",path,"_log2FGene.png"),PlotGenelogFacrossDay2,base_height = 10,base_width =14 )
    
    png(file=paste0(odir,analyse,"PathwayEnrichedIn",numberCond,"DayTypeCell/",path,"_heatmape.png"))
    draw (mH, column_title = path,
          heatmap_legend_side = "bottom", padding = unit(c(2, 2, 2, 60), "mm"))
    dev.off() 
    
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


#####################################################
#### Looking pathways with hierarchical point of view
#####################################################

if (file.exists(tablePathwayHierarchieFull) == F ){
ConvertPathwaysName<-read.csv("data/ReactomePathways.txt",sep = "\t",header = F)
ConvertPathwaysName <- ConvertPathwaysName %>% filter(str_detect(V3,"Mus musculus"))

hierachicalPathway<-read.csv("data/ReactomePathwaysRelation.txt",sep = "\t",header = F)
hierachicalPathway<-hierachicalPathway %>% filter(str_detect(V1,"R-MMU")) %>% mutate("Node1"= str_replace_all(str_to_upper(ConvertPathwaysName[na.omit(match(hierachicalPathway$V1,ConvertPathwaysName$V1)),]$V2), "[\\ /\\-\\:]",'_'))  %>% mutate("Node1"= str_replace_all(Node1,"[,()]", ""))
hierachicalPathway<-hierachicalPathway %>% filter(str_detect(V2,"R-MMU")) %>% mutate("Node2"= str_replace_all(str_to_upper(ConvertPathwaysName[na.omit(match(hierachicalPathway$V2,ConvertPathwaysName$V1)),]$V2),"[\\ /\\-\\:]",'_')) %>% mutate("Node2"= str_replace_all(Node2,"[,()]", ""))

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
  size=ifelse(node2 == "root2", 1,GSEAtempo %>% filter(pathway==paste0("REACTOME_",node2)) %>% select(size) %>% unlist() %>% as.character() %>% as.numeric() )
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

data_edged_zoom2<-data_edged_zoom %>% filter(Node2!="root1") %>% select(Node1,Node2)

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
data_nested[is.na(data_nested)] <- ""
data_nested<-data_nested %>% arrange(desc(level13),desc(level12),desc(level11),desc(level10),desc(level9),desc(level8),desc(level7),desc(level6),desc(level5),desc(level4),desc(level3),desc(level2),desc(level1))

listpathway<-listpathwaytotal

data_nested_print<-data.frame()
data_nested_full<-data.frame()
for(l in rownames(data_nested)){
  newrow=data_nested[l,]
  newrowFull=data_nested[l,]
  pathTF=F
  for (c in 13:1){
    if (data_nested[l,c]!= ""){
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
write.table(data_nested_3, file=paste0( odir,analyse,"ALL_table_REACTOME_Hierachie_print.txt"))

data_nested_full[is.na(data_nested_full)] <- ""
t="sCs"
for (t in Typecellv){
  for (d in daysv[[t]]){
    GSEAtempo<-GSEAsigni%>% filter(type==t,day==d) %>% mutate( pathway2 = str_replace_all(pathway, "REACTOME_", "") )
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
                sizetempo=data_edged_zoom %>% filter(Node2 == data_nested_full[l,c] ) %>% select( size ) %>% unlist() %>% as.character() %>% as.numeric()
                sizetempo=max(sizetempo,na.rm = T )
                }
              pathFound=T 
              newrow[[c]]<- str_replace_all(data_nested_full[l,c],'_',' ')
            }
            else{
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
      
      if (dir.exists(paste0(odir,analyse,d,"_",t,"/")) == F) {
        dir.create(paste0(odir,analyse,d,"_",t,"/"))
      }
      write.table(data_nested_2, file=paste0( odir,analyse,d,"_",t,"/",d,'_',t, "_table_REACTOME_Hierachie.txt"))

      
      data_nested_2$pathString <- paste("roots", data_nested_2$level1, data_nested_2$level2, data_nested_2$level3, data_nested_2$level4,data_nested_2$level5,data_nested_2$level6,data_nested_2$level7, data_nested_2$level8, data_nested_2$level9, data_nested_2$level10,data_nested_2$level11,data_nested_2$level12, sep = "/")
      
      
      data_Node <-as.Node(data_nested_2)
      if (t ==""){
        p <- circlepackeR(data_Node, size = "value",color_min = paste0("hsl(56,80%,80%)"), color_max = paste0("hsl(341,30%,30%)"), width = 1800, height = 1000)
        
      }else{
      p <- circlepackeR(data_Node, size = "value",color_min = paste0("hsl(",hsl_color[[t]],",80%,80%)"), color_max = paste0("hsl(",hsl_color[[t]],",30%,30%)"), width = 1500, height = 1000)
      }
      p
      
      library(htmlwidgets)
      saveWidget(p, file=paste0( odir,analyse,d,"_",t,"/",d,'_',t, "_REACTOME_Hierachie.html"))
      
    }      
  }
}

for (t in Typecellv){
  for(d in daysv[[t]]){
    GSEAtempo<-GSEAsigni%>% filter(day==d,type==t) %>% mutate( pathway2 = str_replace_all(pathway, "REACTOME_", "") )
    if (dim(GSEAtempo)[1] >0){
      rmarkdown::render(paste0('/home/bioinfo/BulkAnalysis_plusNetwork2/scripts/FigureSupp_gsea_Hierarchie.Rmd'),params = list("title"=paste0("Reactome_hierarchie_in_",t,"_",d),"d"=d, "t"=t,"GSEAsigni"=GSEAtempo), output_file = paste0(odir,analyse,d,"_",t,"/",d,"_",t,"_Hierarchical_Figsupp") )
    }
  }
}
for (t in Typecellv){
    GSEAtempo<-GSEAsigni%>% filter(type==t) %>% mutate( pathway2 = str_replace_all(pathway, "REACTOME_", "") )
    if (dim(GSEAtempo)[1] >0){
      d=''
      rmarkdown::render(paste0('/home/bioinfo/BulkAnalysis_plusNetwork2/scripts/FigureSupp_gsea_Hierarchie.Rmd'),params = list("title"=paste0("Reactome_hierarchie_in_",t,"_",d),"d"=d, "t"=t,"GSEAsigni"=GSEAtempo), output_file = paste0(odir,analyse,d,"_",t,"/",d,"_",t,"_Hierarchical_Figsupp") )
  }
}
t="sCs"
for(d in daysv[[t]]){
    GSEAtempo<-GSEAsigni%>% filter(day==d,type==t) %>% mutate( pathway2 = str_replace_all(pathway, "REACTOME_", "") )
    if (dim(GSEAtempo)[1] >0){
      rmarkdown::render(paste0('/home/bioinfo/BulkAnalysis_plusNetwork2/scripts/FigureSupp_gsea_Hierarchie.Rmd'),params = list("title"=paste0("Reactome_hierarchie_in_",t,"_",d),"d"=d, "t"=t,"GSEAsigni"=GSEAtempo), output_file = paste0(odir,analyse,d,"_",t,"/",d,"_",t,"_Hierarchical_Figsupp") )
    }
  }



GSEAtempo<-GSEAsigni%>% filter(day==d,type==t) %>% mutate( pathway2 = str_replace_all(pathway, "REACTOME_", "") )
listpathway<- c(unique(GSEAtempo$pathway2),"NOT_IN_REACTOME_HIERARCHIE","IN_REACTOME_HIERARCHIE")

data_nested_2<-data.frame()
sizeVector<-c()
dico_path<-data.frame()
dico_path<-rbind(dico_path,c("IN_REACTOME_HIERARCHIE",rep(0,13)))
colnames(dico_path)<-c("path",1:13)
pathwayFound<-c()

for(l in rownames(data_nested)){
  newrow=rep("", 13)
  sizetempo=0
  pathFound=F
  for (c in 13:1){
    if (data_nested[l,c]!= ""){
      if ( data_nested[l,c] %in% listpathway){
        if (pathFound == F){
          sizetempo=data_edged_zoom %>% filter(Node2 == data_nested[l,c] ) %>% select( size ) %>% unlist() %>% as.character() %>% as.numeric()
        }
        pathFound=T 
        newrow[[c]]<- str_replace_all(data_nested[l,c],'_',' ')
      }
      else{
        if (pathFound == T){
          if (data_nested[l,c] %in% dico_path$path){
            newrow[[c]]<- paste0(c,'_',dico_path[dico_path$path == data_nested[l,c], c+1])
          }
          else {
            newKey<-dico_path[length(dico_path$path),]
            newKey[["path"]]<-data_nested[l,c]
            newKey[[c+1]]<-as.numeric(newKey[[c+1]])+1
            newrow[[c]]<- paste0(c,'_',newKey[[c+1]])
            dico_path<-rbind(dico_path,newKey)
          }
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
data_edged=data.frame()
    data_edged<-explorer2(hierachicalPathway2,"root2","root2",listSommetparcouru,data_edged,0)
    colnames(data_edged)<-c("Node1","Node2","size","NamePath","level")
    Transfosize=c(max(as.numeric(data_edged$level)):1)
    names(Transfosize)<-c(1:max(as.numeric(data_edged$level)))
    print(Transfosize)
    data_edged <- data_edged %>% mutate ( size = Transfosize[data_edged$level])
    edges <- data_edged %>% select(Node1,Node2)
    vertices <- data_edged%>% mutate(showlabel=ifelse(Node2 %in% Node1 , NamePath,NA)) %>% select (Node2,level,NamePath,size,showlabel) %>% rbind(c("root2",1,"root",7,NA))
    # Rebuild the graph object
    
    mygraph <- graph_from_data_frame( edges, vertices=vertices )
    
    # left
    circles<-ggraph(mygraph, layout = 'circlepack', weight=as.numeric(size) ) + 
      geom_node_circle(aes(fill = level)) +
      geom_node_label(aes(label=showlabel, size=as.numeric(size)),repel = T) +
      theme_void() + 
      theme(legend.position=c(fill="left")) + 
      scale_fill_viridis(discrete = T)
    print(circle)






###
innerplot_fun <- function(fgsea.dayhere, d, CT, NumP, msigdbr_list){
  print("running innerplot_fun")
  fgseaRes <- fgsea.dayhere[[CT]] %>% mutate(log2err=replace_na(log2err, 1))  ## not use drop_na : all down faps D2 are log2err NA
  print(min(fgseaRes$padj))
  topPathwaysUp <- fgseaRes[ES>0][head(order(padj),n=NumP), ]
  print(paste("nb of upregul paths", length(topPathwaysUp)))
  topPathwaysDown <- fgseaRes[ES<0][head(order(padj), n=NumP), ]
  print(paste("nb of downreg paths", length(topPathwaysDown)))
  topPathBoth <- c(topPathwaysUp$pathway, rev(topPathwaysDown$pathway))
  a <- sort(getgeneslist_mod(topPathwaysUp, d), decreasing = T)
  b <- sort(getgeneslist_mod(topPathwaysDown, d), decreasing = T)
  gseagenes <- c(a, b)
  gseagenes <- sort(gseagenes, decreasing=T)
  print(gseagenes)
  ouif = fgsea::plotGseaTable(msigdbr_list[topPathBoth], gseagenes, fgseaRes, 
                              gseaParam = 0.5 ,  colwidths = c(3, 3, 0.8, 1.2, 1.2),  render=F) 
  print('doneplotgseatable')
  plotsenrichu_ <- list()
  for (i in 1:NumP){
    pup = topPathwaysUp[i,]
    print("a path up")
    print(pup)
    tmpup <- NULL
    tryCatch({ 
      tmpup <- fgsea::plotEnrichment(msigdbr_list[[pup$pathway]], gseagenes) + 
        labs(title = str_replace(pup$pathway, "REACTOME_", ""),
             caption=paste( "(NES:", round(pup$NES, 2), ", padj :", 
                            round(pup$padj,2),")") )
    }, 
    warning = function(w) {}, error = function(e){print("paths up error")} ) # endtrycatch
    if ((dim(tmpup$data)[1] <= 2) || is.null(tmpup) ){
      plotsenrichu_[[i]] <- NULL
    }else{ plotsenrichu_[[i]] <- tmpup + theme(plot.title = element_text(size=8))}
  }
  print("done enrichu")
  plotsenrichdw_ <- list()
  for (i in 1:NumP){
    pdw = topPathwaysDown[i,]
    print("b path down")
    print(pdw)
    tmpdw <- NULL
    tryCatch({
      tmpdw <- fgsea::plotEnrichment(msigdbr_list[[pdw$pathway]], stats=gseagenes) + 
        labs(title = str_replace(pdw$pathway,"REACTOME_", ""),
             caption=paste( "(NES:", round(pdw$NES,2), ", padj :", 
                            round(pdw$padj,2),")") )
    }, 
    warning = function(w){}, error = function(e){print("paths down error")} ) # trycatch
    if (dim(tmpdw$data)[1] <= 2 || is.null(tmpdw)){
      plotsenrichdw_[[i]] <- NULL
    } else {plotsenrichdw_[[i]] <- tmpdw + theme(plot.title = element_text(size=8)) }
  }
  print("done down")
  print(
    plot_grid(
      # inner plotgrid 1 : the 'plotGseaTable',
      plot_grid(ggdraw() + draw_label(paste(d, CT, ": Top enriched Pathways (GSEA), Old vs Young")),
                plot_grid(NULL, ouif, NULL, nrow=1, rel_widths = c(10,9,1)), # NULL elem helps de-truncate names
                nrow=2, rel_heights = c(2, 10)) ,
      # inner plotgrid 2 : the 'plotEnrichment' : u_ == UP, dw_ == DOWN
      plot_grid(
        plot_grid(plotlist = plotsenrichu_, nrow = 1, ncol=NumP+1),
        plot_grid(plotlist = plotsenrichdw_,  nrow = 1, ncol=NumP+1) , nrow = 2, rel_heights = c(5,5)), 
      
      nrow=2, rel_heights = c(4,5)
    )# end plotgrid all  
  ) # end print
  print("done plotgrid")
}

getgeneslist_mod <- function(gseadatafr, d){
  outi <- c()
  for (i in gseadatafr$leadingEdge){
    outi <- c(outi, i)
  }
  moo <- unique(outi)
  tmpdfdeg <- fullDEsta %>% filter(day==d)
  lfcs <- tmpdfdeg[match(moo, tmpdfdeg$symbol),]$log2FoldChange
  names(lfcs) <- moo
  return(lfcs[!is.na(lfcs)])
}

plotme_mod <- function(pathsFiltered, d, outfileprefix, NumP=5){
  fgsea.dayhere <- pathsFiltered[[d]]
  print(paste("pathways object is a list of lists, inner list has celltypes:",
              ifelse((typeof(fgsea.dayhere) == "list" & 
                        is.null(colnames(fgsea.dayhere))), "ok",  #colnames(a list) is null
                     "ERROR no celltypes in names(fgsea.dayhere), type prblm" )))
  types = names(fgsea.dayhere)
  print(types)
  pdf(paste0(odir, "GSEA/", outfileprefix, d,".pdf"), width = 14, height=8)
  for (CT in types){
    print(CT)
    print( innerplot_fun( fgsea.dayhere, d, CT, NumP, msigdbr_list) )
  }  
  dev.off()
}
rm(d)
plotme_mod(pathsFiltered, 'D0', "result_", NumP=5)
plotme_mod(pathsFiltered, 'D2', "result_", NumP=5)
plotme_mod(pathsFiltered, 'D4', "result_", NumP=5) 
plotme_mod(pathsFiltered, 'D7', "result_", NumP=5) 



# ============= Prepare matrices for individual (dummy) plots ==================
## prepare matrices for heatmaps (small heatmaps )
mat4heatmapneeded<-TRUE
if (mat4heatmapneeded){
  m_ = list()
  for (k in c('D0','D2', 'D4', 'D7')){
    cts <- unique(names(pathsFiltered[[k]])) 
    m_[[k]] <- list()
    for (CT in cts){
      m_[[k]][[CT]] <- list("UP"=NULL,"DOWN"=NULL)
      for(S in c("UP","DOWN")){
        here.data = as.data.frame(pathsFiltered[[k]][[CT]])
        here.data <- here.data %>% filter(sens == S) %>% 
          mutate(path4graph = str_replace(pathway, "REACTOME_", "")) 
        here.dico = list()
        if(length(here.data$path4graph) == length(unique(here.data$path4graph) )){
          print(paste("ok, processing", k, CT, S))
          rownames(here.data) =  here.data$path4graph
          aggreg_genes = c()
          for (path in here.data$path4graph){
            tmpstr <- str_replace_all(here.data[path,]$leadingEdge, 
                                      c('c\\("' = '',   '"' = '', '\\)' = ''))
            here.dico[[path]] <- unname(unlist(str_split(tmpstr,", ")))
            aggreg_genes <- c(aggreg_genes, here.dico[[path]] )
          }
          colnamesgenes = unique(aggreg_genes)
          matrx = array(NA, dim=c(length(here.data$path4graph), length(colnamesgenes)))
          rownames(matrx) = here.data$path4graph
          colnames(matrx) = colnamesgenes
          for (path in here.data$path4graph){
            localgenes <- here.dico[[path]]
            for (g in localgenes){
              matrx[path,g] <-  here.data[path,]$NES # if more than once, yields 2 or more
            }
          }
          # clear rows and columns containing 2 or less values
          matrx <- matrx[rownames(matrx) != "DISEASE",] # take away "DISEASE" ! 
          rtokeep <- apply(matrx, 1, function(x) sum(!is.na(x)) >= 3)
          newm <- matrx[rtokeep,]
          m_[[k]][[CT]][[S]] <- newm
        } else {
          print(paste("error, this", k, CT ,"  contains repeated pathwaynames") )
          stop()
        } # endif pathways vectors are unique (no elements repeated)
      } # end for sens UP or DOWN
    }# end for CT in cts
  }# end for k in vector of days
  saveRDS(m_, file=paste0(odir,"GSEA/rds/fgsea_matrices4heatmaps.rds" ))
} # end if mat4heatmapneeded

# ===================== Prepare giant single matrix ============================
#           and exclude DISEASE (use str_detect) and DEVELOPMENT 
nbPathwaysDisease<-0
givemeltedDataframe <- function(pathsFiltered, minnbg, nbminpadj, COLLECTION){
  r_path = c(); r_sens = c(); r_NES = c(); 
  r_gene = c(); r_celltype = c(); r_day = c(); r_padj = c()
  counter = 0
  for (d in names(pathsFiltered)){
    for (t in names(pathsFiltered[[d]])){
      tmp <- pathsFiltered[[d]][[t]] %>% filter(str_detect(pathway, COLLECTION)) %>%
        group_by(sens) %>% slice_min(padj, n=nbminpadj)
      path_ = tmp$pathway
      for (i in 1:length(path_)){
        pw = tmp[i,]$pathway
        NES = tmp[i,]$NES
        se = tmp[i,]$sens
        pj = tmp[i,]$padj
        if (str_detect(pw, "DISEASE") | str_detect(pw, "DEVELOPMENT") | str_detect(pw, "INFECTION") ){
          nbPathwaysDisease<-nbPathwaysDisease+1
          print(paste("excluded disease, devel and infection terms",nbPathwaysDisease))
        } else {
          for(lg in tmp[i,]$leadingEdge){
            if(length(lg) >= minnbg ){
              print(lg)
              counter = counter + 1
              for (g in lg){
                r_path = c(r_path, pw)
                r_sens = c(r_sens, se) 
                r_NES = c(r_NES, NES)
                r_gene =  c(r_gene, g)
                r_celltype = c(r_celltype, t)
                r_day = c(r_day, d)
                r_padj = c(r_padj, pj)
                print(counter)
              } 
            }
          }
        }
      }#end for i in 1:length(path_)
    }
  }
  dfresu <- data.frame("pathway" = r_path, "sens"= r_sens, "NES" = r_NES,
                       "gene" = r_gene, "celltype" = r_celltype,
                       "day"=r_day, "padj" = r_padj)
  dfresu <- dfresu %>% mutate(unigene = paste0(gene,"_", celltype),
                              unipath = paste0(pathway,"_", day)) 
  
  return(dfresu)
} 

MELTED <- givemeltedDataframe(pathsFiltered, 5, 15, "REACTOME") 
names(MELTED) 
head(MELTED)
# path sens        NES   gene celltype day    unigene                                            unipath
# 1 REACTOME_ABC_FAMILY_PROTEINS_MEDIATED_TRANSPORT   UP  0.7683668  Abcb4      ECs  D0  Abcb4_ECs REACTOME_ABC_FAMILY_PROTEINS_MEDIATED_TRANSPORT_D0
# 2 REACTOME_ABC_FAMILY_PROTEINS_MEDIATED_TRANSPORT   UP  0.7683668  Abcc1      ECs  D0  Abcc1_ECs REACTOME_ABC_FAMILY_PROTEINS_MEDIATED_TRANSPORT_D0
# 3 REACTOME_ABC_FAMILY_PROTEINS_MEDIATED_TRANSPORT   UP  0.7683668  Abcc9      ECs  D0  Abcc9_ECs REACTOME_ABC_FAMILY_PROTEINS_MEDIATED_TRANSPORT_D0
# 4 REACTOME_ABC_FAMILY_PROTEINS_MEDIATED_TRANSPORT   UP  0.7683668  Abca5      ECs  D0  Abca5_ECs REACTOME_ABC_FAMILY_PROTEINS_MEDIATED_TRANSPORT_D0
# 5                 REACTOME_ADAPTIVE_IMMUNE_SYSTEM DOWN -1.0629117 Zbtb16      ECs  D0 Zbtb16_ECs                 REACTOME_ADAPTIVE_IMMUNE_SYSTEM_D0

print("adding information about genes (padj, lfc), from 'rds/shot_rds_full.rds'")
toto <- MELTED %>% select(unigene, gene, celltype, day)
fuu <- fullDEsta %>% select(symbol, padj, log2FoldChange, type, day) %>% 
  mutate(gene = symbol, celltype = type, genepadj = padj)
tototo <- left_join(toto, fuu, by = c("day", "gene", "celltype" ) )
tototo <- tototo  %>% select(unigene, genepadj, log2FoldChange, day )
ttp = left_join(MELTED, tototo, by=c("unigene", "day"))  %>% unique() 
MELTED = ttp

rm(toto, fuu, tototo, ttp)

approve <- MELTED %>% group_by(celltype, day) %>% slice_min(genepadj, n=10) %>%
  filter(abs(log2FoldChange)>=0.3)
length(unique(approve$pathway))
length(unique(approve$unigene))

MELTED <- approve
# funciton to build megamatrix:
dfTomatrix <- function(dfx){
  mat <- array(NA,dim=c( length(unique(dfx$unipath)),
                         length(unique(dfx$unigene))
  ))
  rownames(mat) <- unique(dfx$unipath)
  colnames(mat) <- unique(dfx$unigene)
  for (p in dfx$unipath){
    for (g in dfx$unigene){
      tryCatch({
        val <- dfx[dfx$unipath == p & dfx$unigene == g,]$NES
      },
      warning = function(w){
      },
      error = function(e){
      })
      if (length(val)!= 0 ){
        mat[p,g] <- val
        print(val)
      }
    }
  }
  return(mat)
}

print("building megamatrix, take a cup of coffee (takes 15 minutes)")
m_reac <- dfTomatrix(MELTED)
print("clearing megamatrix")
keepro <- apply(m_reac, 1, function(r) sum(!is.na(r))>1)
m_reac2 <- m_reac[keepro,]

reac_daysvec <- sapply(rownames(m_reac2), function(x) {
  ltmp <- unlist(str_split(x, "_")) 
  return(ltmp[length(ltmp)])})
reac_celltyvec <- sapply(colnames(m_reac2),function(x){
  ltmp <- unlist(str_split(x, "_"))
  return(ltmp[length(ltmp)])
})
reac_newrownames <-  make.unique( sapply( names(reac_daysvec), function(x){  
  # "REACTOME_COLLAGEN_FORMATION_D7"
  ltmp <- unlist(str_split(str_replace(x,"REACTOME_",""), "_"))  
  mot <- paste0(ltmp[1:(length(ltmp)-1)], collapse="_")
  return(mot)  # "COLLAGEN_FORMATION
}) )
reac_newcolnames <- make.unique( sapply(names(reac_celltyvec), function(x){
  ltmp <- unlist(str_split(x, "_"))  # "Lamc1_sCs"
  return(ltmp[1]) # Lamc1
}) )
# assigne new names to all objects
rownames(m_reac2) <- reac_newrownames
colnames(m_reac2) <- reac_newcolnames
names(reac_daysvec) <- reac_newrownames
names(reac_celltyvec) <- reac_newcolnames

reac_splitdays <- data.frame(x=reac_daysvec)
rownames(reac_splitdays) <- rownames(m_reac2)
reac_splitcols <- data.frame(y=reac_celltyvec)
rownames(reac_splitcols) <- colnames(m_reac2)

ha = HeatmapAnnotation(foo = anno_block(gp = gpar(fill = c("darkblue", "orange",
                                                           "firebrick","violet",
                                                           "darkgreen","royalblue") )
)) #  labels = c("ECs","FAPs","M1","M2","Neutro","sCs") 
oh <- ComplexHeatmap::Heatmap(m_reac2, 
                              na_col = "whitesmoke",
                              cluster_rows = FALSE,
                              cluster_columns = FALSE,
                              row_split = reac_splitdays,
                              column_split = reac_splitcols,
                              row_names_gp = gpar(fontsize = 9),
                              column_names_gp = gpar(fontsize = 8),
                              name = "REACTOME",
                              column_names_rot = 45,
                              heatmap_legend_param = list(title="NES",
                                                          direction = "horizontal"),
                              top_annotation = ha,
                              use_raster = TRUE
)

pdf(paste0(odir,"GSEA/fgsea_byday_bycelltype.pdf"), width = 18, height = 10)
reactomeGSEAheatmap<-draw (oh, column_title = "GSEA Old vs Young (REACTOME)",
                           heatmap_legend_side = "bottom", padding = unit(c(2, 2, 2, 60), "mm"))
reactomeGSEAheatmap
save_plot(paste0(odir,"fgsea_byday_bycelltype.png"),reactomeGSEAheatmap)
dev.off()
