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
library(ggplot2)
setwd("~/BulkAnalysis_plusNetwork/")
analyse<-"analysis/AnalysisGSEA/"
odir <- "exam_INTER_conditions/static/"
gseardsFull = "GSEAfull/rds/"
tablesneeded <- F
mat4heatmapneeded <- F
gseaoutfull = "fgsea_full.rds" 
gseaoutfiltered = "fgsea_filtered.rds"

fullrds <- readRDS(paste0(odir,"GSEAfull/rds/", gseaoutfull))
pathsFiltered = readRDS( paste0(odir, "GSEAfull/rds/",  gseaoutfiltered) )

thegmt <- msigdbr(species = "Mus musculus", 
                  category = 'C2', 
                  subcategory=c('CP:REACTOME'))
msigdbr_list = split(x = thegmt$gene_symbol, f = thegmt$gs_name)

fullDEsta = readRDS(paste0(odir, "rds/shot_rds_full.rds"))
daysv<-unique(fullDEsta$day)
Typecellv<-unique(fullDEsta$type)
####Extract individual table #####
if (tablesneeded){
  system(paste0("cd ",odir,"GSEAfull/; ", 
                "if [ ! -d csv/ ]; then mkdir csv; else echo 'csv/ exists, no action';fi"))
  print(" ** see keys in full results rds list of lists ** ")
  print(names(fullrds)) # the days
  print(names(fullrds[["D0"]])) # the cell types (that exist in given day)
  days = names(fullrds)
  OUT <- createWorkbook()
  for (d in names(fullrds)){
    print(paste("\n^^^^^", d, "^^^^^"))
    for (t in names(fullrds[[d]])){
      print(paste("saving", d, t, "as csv"))
      print(t)
      print(dim(fullrds[[d]][[t]]))
      tmpdf <- fullrds[[d]][[t]]
      
      tmpdf <- tmpdf %>% 
        mutate(gene_symbols = sapply(leadingEdge, function(x) paste(x, collapse=", ")),
               .before = leadingEdge) %>%
        select(-leadingEdge)
      tmpdf <- tmpdf %>% arrange(padj,desc(abs(NES)),desc(length(gene_symbols)))
      
      write.table(tmpdf, paste0(odir,"GSEAfull/csv/", t,"_",d,"_pathways.csv"), sep="\t", 
                  col.names = T, row.names=F              )
      addWorksheet(OUT, paste0(t,'_',d))
      writeData(OUT,sheet= paste0(t,'_',d),x=tmpdf)
    }
  }
  saveWorkbook(OUT,paste0(odir,"GSEAfull/csv/Total_GSE_pathways.xlsx"),overwrite = T)
  rm(fullrds)
}
####Function color
darker.col = function(color, how.much = 30){
  colorRampPalette(c("white",color,"black"))(9)[how.much]
}
tablecolor<-data.frame("ECs"=darker.col("#0072B2", 4:7),"FAPs"=darker.col("#F0E442", 3:6),"M1"=darker.col("#D55E00", 3:6),"M2"=darker.col("#CC79A7", 3:6),"Neutro"=darker.col("#009E73", 3:6),"sCs"=darker.col("#56B4E9",2:5))
rownames(tablecolor)<-c("D0","D2","D4","D7")

### 
# Extract summary signi pathway enrichie
###

XX_full = readRDS(paste0(odir,gseardsFull, gseaoutfull))
padjofpathsNumbersigni_full = array(NA, dim=c(4,6))
padjofpathsDEsigni = array(NA, dim=c(4,6))

rownames(padjofpathsNumbersigni_full) =  names(XX_full)
colnames(padjofpathsNumbersigni_full) = names(XX_full[['D2']])

rownames(padjofpathsDEsigni) = names(XX_full)
colnames(padjofpathsDEsigni) = names(XX_full[['D2']])

for (d in names(XX_full)){
  for (n in names(XX_full[[d]])){
    print(length(unique(XX_full[[d]][[n]][XX_full[[d]][[n]]$padj <= 0.05,]$pathway)))
    padjofpathsNumbersigni_full[d, n] <-length(unique(XX_full[[d]][[n]][XX_full[[d]][[n]]$padj <= 0.05,]$pathway))
    padjofpathsDEsigni[d, n] <- length(fullDEsta %>% filter(day == d & type == n & padj <=0.05) %>% dplyr::select(symbol)  %>% unlist() %>% unique())
  }
}
totalDay = padjofpathsNumbersigni_full %>% rowSums(na.rm = TRUE)
toltalTypeCell = padjofpathsNumbersigni_full %>% colSums(na.rm = TRUE)
totalDay = c(totalDay,sum(totalDay))
padjofpathsNumbersigni_full = rbind(padjofpathsNumbersigni_full, "TotalPathwaySigni"= toltalTypeCell)
padjofpathsNumbersigni_full = cbind(padjofpathsNumbersigni_full, "TotalPathwaySigni"= totalDay)
padjofpathsNumbersigni_full
####
# Table GSEA concat to explore combinaison pathway
####
GSEAconcat<-tibble() # table concat to explore combinaison pathway signi
for (d in daysv){
  for (n in names(XX_full[[d]])){
    GSEAconcat<-rbind(GSEAconcat, XX_full[[d]][[n]] %>% mutate(day=d) %>% mutate(type=n))
  }
}
GSEAsigni<- GSEAconcat %>% filter(padj<=0.05) %>% arrange(day,type) %>% mutate(sens=ifelse(NES < 0,"down", "up"))
RepartitionPathwaySummary<-table(table(GSEAsigni$pathway))
RepartitionPathway<-as.data.frame(table(GSEAsigni$pathway)) 
library(oce)
library(ggraph)
####Function color
darker.col = function(color, how.much = 30){
       colorRampPalette(c("white",color,"black"))(9)[how.much]
   }
tablecolor<-data.frame("ECs"=darker.col("#0072B2", 4:7),"FAPs"=darker.col("#F0E442", 3:6),"M1"=darker.col("#D55E00", 3:6),"M2"=darker.col("#CC79A7", 3:6),"Neutro"=darker.col("#009E73", 3:6),"sCs"=darker.col("#56B4E9",2:5))
rownames(tablecolor)<-c("D0","D2","D4","D7")

daysv<-list()
daysv$ECs=c("D0","D2","D4","D7")
daysv$FAPs<-c("D0","D2","D4","D7")
daysv$M1<-c("D2","D4")
daysv$M2<-c("D2","D4","D7")
daysv$Neutro<-c("D2")
daysv$sCs<-c("D0","D2","D4","D7")
Typecellv<-vectTypeCell <- c("ECs","FAPs","M1","M2","Neutro","sCs")

DynamicPathwayAcrossCond <- function(RepartitionPathway,GSEAsigni,n){
  PathwaySigniCond<-RepartitionPathway %>% filter(Freq==n)
  PathwaySigniCondtable<- lapply(PathwaySigniCond$Var1 , function(x) dplyr::filter(GSEAsigni, pathway == x))
  DynamicCond <-lapply(PathwaySigniCondtable, function(x) paste0(x$day,'_',x$type,'_',x$sens) %>% sort())
  names(DynamicCond)<-lapply(PathwaySigniCondtable, function(x) unique(x$pathway) )
  DynamicCondconcat<-c()
  edgelist<-data.frame()
  for(Dym in DynamicCond){DynamicCondconcat=c(DynamicCondconcat,str_c(Dym,collapse = "."))
  }
  names(DynamicCondconcat)<-lapply(PathwaySigniCondtable, function(x) unique(x$pathway))
  DynamicCombi =table(DynamicCondconcat)
  return(DynamicCondconcat)
}
      

GetPlot<-function(Dynamic6Condconcat,GSEAconcat, nbCond ,vertical){
  tableDynamique2<-sort(table(Dynamic6Condconcat),decreasing = T)
  tablePlot<-data.frame()
  for ( name in names(tableDynamique2)){
    edgeslabel<-unlist(str_split(name, pattern = "\\."))
    for ( cond in edgeslabel){
      infoCond<-unlist(str_split(cond, pattern = "_"))
      print(infoCond)
      print(infoCond[1])
      dayt<-infoCond[1]
      typet<-infoCond[2]
      pathways<-as.tibble(Dynamic6Condconcat) %>% mutate(pathway = names(Dynamic6Condconcat)) %>% filter(value == name )%>% select(pathway) %>% unlist()
      numberPath<-length(pathways)
      for (path in 1:length(pathways) ){
        print(pathways[[path]])
        padj<-GSEAconcat %>% filter(pathway == pathways[[path]]) %>% filter(day == dayt) %>% filter(type == typet) %>% select(padj) %>% unlist() %>% as.character()
        NES<-GSEAconcat %>% filter(pathway == pathways[[path]]) %>% filter(day == dayt) %>% filter(type == typet) %>% select(NES) %>% unlist() %>% as.character()
        name2=paste0(str_replace_all(name, "\\.", "  "),' ,there are ',numberPath,' with this dynamic')
        tablePlot<-rbind(tablePlot,c(name2,pathways[[path]],dayt,typet,padj,NES,numberPath,tablecolor[dayt,typet]))
      }
    }
  }
  names(tablePlot)<-c("name","pathway","day","type","padj","NES","numberPath","color")
  tablePlot <- tablePlot %>% mutate( condition = paste0(day,'_',type))
  col<- tablePlot %>% arrange(condition) %>% select(color) %>% unlist() %>% as.character() %>% unique()
  tablePlot <-tablePlot %>% mutate( name = as.factor(name)) %>% mutate(name = fct_reorder(name, as.numeric(numberPath), .desc = T) )
  typebytime<-ggplot(tablePlot, aes(x=day,y=as.numeric(NES),fill=condition,color=condition,group=paste0(pathway,'_',type))) +
    geom_point(aes(size=5,alpha=0.8),show.legend = FALSE)+
    geom_line()+
    scale_color_manual(values=col)+
    ggtitle( paste0("Pathways with ",nbCond," conditions enriched")) +
    facet_wrap(~name,dir = "v",ncol=vertical) +
    theme_minimal()
  print(typebytime)
  return(tablePlot)
}


Dynamic6Condconcat=DynamicPathwayAcrossCond(RepartitionPathway,GSEAsigni,6)
Dynamic5Condconcat=DynamicPathwayAcrossCond(RepartitionPathway,GSEAsigni,5)
Dynamic4Condconcat=DynamicPathwayAcrossCond(RepartitionPathway,GSEAsigni,4)
Dynamic3Condconcat=DynamicPathwayAcrossCond(RepartitionPathway,GSEAsigni,3)
Dynamic2Condconcat=DynamicPathwayAcrossCond(RepartitionPathway,GSEAsigni,2)
GetPlot(Dynamic6Condconcat,GSEAconcat,6, 1)
GetPlot(Dynamic5Condconcat,GSEAconcat,5,1)
GetPlot(Dynamic4Condconcat,GSEAconcat,4,1)
GetPlot(Dynamic3Condconcat,GSEAconcat,3, NA)
GetPlot(Dynamic2Condconcat,GSEAconcat,2, NA)
