Organize GSEA results into tables
#   check fgsea_dayandcelltype.R
# output dir: exam_INTER_conditions/static/GSEA/ 
# and pdf plots
# --
# johaGL
library(tidyverse)
library(ComplexHeatmap)
library(msigdbr)
library(cowplot)


#Directories
setwd("~/BulkAnalysis_plusNetwork/")
dir<-getwd()
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



####################
####Function colors
####################

orderTypecell=c("ECs","FAPs","MuSCs","Neutrophils","Inflammatory-Mac","Resolving-Mac")
colorsType=c("#10b387ff","#3d85c6ff","#b171f1ff","#f0e442ff","#ff9900ff","#cc0000ff")
names(colorsType)=orderTypecell



if (mat4heatmapneeded){
  m_ = list()
  for (k in c('D0','D2', 'D4', 'D7')){
    cts <- unique(names(filteredGSEA[[k]])) 
    m_[[k]] <- list()
    for (CT in cts){
      m_[[k]][[CT]] <- list("UP"=NULL,"DOWN"=NULL)
      for(S in c("UP","DOWN")){
        here.data = as.data.frame(filteredGSEA[[k]][[CT]])
        here.data <- here.data %>% dplyr::filter(sens == S) %>% 
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
          matrx <- matrx[rownames(matrx) != "DISEASE",] # take away "DISEASE" ! 
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
} 