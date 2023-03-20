####################
# Analysis NatmiOut
#####################

library(tidyverse)
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
library(ggplot2)
library(ggExtra)
library(arcdiagram)
library(plyr)
library(circlepackeR) 
library(igraph)
library(data.tree)
library(spaceNtime)
library("readxl")

##################################
#### Path files to load or created
##################################

#Directories
setwd("~/BulkAnalysis_plusNetwork/")
analyseDiff_TPM<-"natmiOut_TPM/analysis/Diff/"
analyseOld_TPM<-"natmiOut_TPM/analysis/Old/"
analyseYoung_TPM<-"natmiOut_TPM/analysis/Young/"
natmiOut_TPM<-'natmiOut_TPM/'
analyseDiff_CN<-"natmiOut_CountNormalised/analysis/Diff/"
analyseOld_CN<-"natmiOut_CountNormalised/analysis/Old/"
analyseYoung_CN<-"natmiOut_CountNormalised/analysis/Young/"
natmiOut_CN<-'natmiOut_CountNormalised/'
EdgesNatmi<-"Network_exp_0_spe_0_det_0.6_top_0_signal_lrc2p_weight_mean/"
DeltaNatmi<-"Delta_edges_lrc2p/"
odir <- "/home/bioinfo/BulkAnalysis_plusNetwork/"
dataGene = "exam_INTER_conditions/static/rds/"
dataGSEA="exam_INTER_conditions/static/GSEA/rds/"
HierarchieData= "exam_INTER_conditions/static/GSEA/HierarchieData/"
inDataNatmi<-"inDataNatmi/"

#Lo0ads
DEGfile =  "TableDynamicUpDownDEG.rds"
GSEAfile = "TableGSEAsigniDynamics.rds"
litteratureSupportRL<- "../NATMI/41467_2020_18873_MOESM4_ESM.xlsx"

#Created
concatYoung_TPM="TableRLYoung"
concatOld_TPM="TableRLOld"
concatDiff_TPM="TableDiff"

concatYoung_CN="TableRLYoung"
concatOld_CN="TableRLOld"
concatDiff_CN="TableDiff"

concatExprAllReplicateYoung_TPM="ConcatAllReplicate_TPM_Young.csv"
concatExprAllReplicateYoung_CN="ConcatAllReplicate_CountNormalised_Young.csv"
concatExprAllReplicateOld_TPM="ConcatAllReplicate_TPM_Old.csv"
concatExprAllReplicateOld_CN="ConcatAllReplicate_CountNormalised_Old.csv"

ExprMeanTPMYoung_SumDay="ExpressionMeanTPM_SumDay_Young.csv"
ExpreMeanCNYoung_SumDay="ExpressionMeanCountNormalised_SumDay_Young.csv"
ExprMeanTPMOld_SumDay="ExpressionMeanTPM_SumDay_Old.csv"
ExpreMeanCNOld_SumDay="ExpressionMeanCountNormalised_SumDay_Old.csv"

##################################
#### Load files
##################################

metadata = readRDS("data/metadata.rds")
design  = read.csv("data/design.csv")
DEGinfo =  readRDS(paste0(dataGene,DEGfile) )
fullDEsta = readRDS(paste0(dataGene, "shot_rds_full.rds"))

GSEAinfo = readRDS(paste0(dataGSEA,GSEAfile) )
GSEADynamics<-readRDS(paste0(odir,dataGSEA,"/TableGSEAsigniDynamics.rds"))

LigandLocation<-read_excel(litteratureSupportRL, sheet = "literature_support")
correspondance_9606_10090<-read_csv("/home/bioinfo/NATMI/homology/correspondance_9606_10090.csv" )
LigandLocation_mmu<- merge(LigandLocation,correspondance_9606_10090,by.x="Ligand gene symbol",by.y="9606") 
LigandLocation_mmu<- LigandLocation_mmu %>% select(`10090`,`Ligand location`)


thegmt <- msigdbr(species = "Mus musculus", 
                  category = 'C2', 
                  subcategory=c('CP:REACTOME'))
msigdbr_list = split(x = thegmt$gene_symbol, f = thegmt$gs_name)
genes_df <- read.table("data/genesinfo.csv", sep="\t", header=T)


daysv<-list()
daysv<-lapply( unique(DEGinfo$type), function(t) unique(DEGinfo[DEGinfo$type == t,]$day) )
names(daysv)<-unique(DEGinfo$type)
Typecellv<-unique(DEGinfo$type)
days<-unique(DEGinfo$day)


####################################
#Recap Nb replicat by conddition
####################################
NbReplicatYoung=array(NA, dim=c(length(days),length(Typecellv)))
NbReplicatOld=array(NA, dim=c(length(days),length(Typecellv)))
rownames(NbReplicatYoung) =  days
colnames(NbReplicatYoung) = Typecellv
rownames(NbReplicatOld) =  days
colnames(NbReplicatOld) = Typecellv
for (t in Typecellv){
  for (d in daysv[[t]]){
    #NbReplicatYoung[d,t]<-dim(metadata %>% filter(time==d,type==t,age=="Young"))[1]
    #NbReplicatOld[d,t]<-dim(metadata %>% filter(time==d,type==t,age=="Old"))[1]
    NbReplicatYoung[d,t]<-dim(design %>% filter(time==d,type==paste0(t,'-RNA'),age=="Young"))[1]
    NbReplicatOld[d,t]<-dim(design %>% filter(time==d,type==paste0(t,'-RNA'),age=="Old"))[1]
  }
}

#### Concat DegInfo & fulllDEsta

GSEADynamics2<-data_frame()
DEGinfo2<-data_frame()

for (t in Typecellv){
  for (d in daysv[[t]]) {
    print(t)
    print(d)
    DEGinfo_dt <-  DEGinfo %>% filter(day==d,type==t)
    GSEADynamics_dt <- GSEADynamics %>% filter(day==d,type==t)
    if (dim(GSEADynamics_dt)[[1]] > 0 ){
      GSEADynamics_dt <- GSEADynamics_dt %>% mutate(leadingEdgeDEG= lapply(GSEADynamics_dt$leadingEdge, function(listg) lapply(listg, function(g) if (g %in% DEGinfo_dt$symbol) { g }) %>% unlist() ))
      GSEADynamics2<-rbind(GSEADynamics2,GSEADynamics_dt)
      
      PathwayGene<-lapply(GSEADynamics_dt$leadingEdge, function(listg) lapply(listg, function(g) if (g %in% DEGinfo_dt$symbol) { g }) %>% unlist() )
      names(PathwayGene)<-GSEADynamics_dt$pathway
      #GenePathways<-lapply(DEGinfo$symbol, function(g) mapply( function(listg,path) if (g %in% listg) {path}, PathwayGene, names(PathwayGene)) %>% unlist() %>% unlist())
      #names(GenePathways)<-DEGinfo$symbol
      tmp = list()
      for ( p in names(PathwayGene)){
        if ( !is.null(PathwayGene[[paste(p)]])){
        tmp[[p]] = data.frame(genes=PathwayGene[[paste(p)]],pathwayname=p)
        }
      }
      tmp2=do.call('rbind',tmp) 
      tmp2<- tmp2 %>% mutate(test="Listpathways")
      GenePathways=tmp2 %>% pivot_wider(names_from = test, values_from = pathwayname) 
      DEGinfo_dt= DEGinfo_dt %>% mutate(listPathways= GenePathways[match(DEGinfo_dt$symbol,GenePathways$genes), ]$Listpathways)
      
      DEGinfo2<-rbind(DEGinfo2,DEGinfo_dt)
    } else{
      DEGinfo_dt<-DEGinfo_dt %>% mutate(listPathways = rep(NA, length(DEGinfo_dt$symbol) ))
      DEGinfo2<-rbind(DEGinfo2,DEGinfo_dt)
    }
  }
}

#Concat NatmiIn to get raw expression of each replicat

if (! file.exists(paste0(inDataNatmi,ExpreMeanCNOld_SumDay)) ) {
  tableExpressionOldTPM<-read.table(file = paste0(inDataNatmi,"TPM_OldD0.txt"),header = T)
  tableExpressionYoungTPM<-read.table(file = paste0(inDataNatmi,"TPM_YoungD0.txt"),header = T)
  for (d in 2:length(days)){
    expresstempo<-read.table(file = paste0(inDataNatmi,"TPM_Old",days[d],".txt"),header = T)
    expresstempo<- expresstempo %>% mutate( "MeanTPM.{d}" := rowSums(.) )
    tableExpressionOldTPM<-merge(tableExpressionOldTPM,expresstempo,by.x="gene",by.y="gene")
    expresstempo<-read.table(file = paste0(inDataNatmi,"TPM_Young",days[d],".txt"),header = T)
    expresstempo<- expresstempo %>% mutate( "MeanTPM.{d}" := rowSums(.) )
    tableExpressionYoungTPM<-merge(tableExpressionYoungTPM,expresstempo,by.x="gene",by.y="gene")
  }
  tableExpressionYoungTPM<- tableExpressionYoungTPM %>% mutate(symbol=genes_df[match(tableExpressionYoungTPM$gene,genes_df$Geneid),]$symbol)
  tableExpressionOldTPM<- tableExpressionOldTPM %>% mutate(symbol=genes_df[match(tableExpressionOldTPM$gene,genes_df$Geneid),]$symbol)
  
  write_csv(tableExpressionYoungTPM,paste0(inDataNatmi,concatExprAllReplicateYoung_TPM))
  write_csv(tableExpressionOldTPM,paste0(inDataNatmi,concatExprAllReplicateOld_TPM))
  
  tableExpressionOldCountNormalised<-read.table(file = paste0(inDataNatmi,"CountNormalised_OldD0.txt"),header = T)
  tableExpressionYoungCountNormalised<-read.table(file = paste0(inDataNatmi,"CountNormalised_YoungD0.txt"),header = T)
  for (d in 2:length(days)){
    expresstempo<-read.table(file = paste0(inDataNatmi,"CountNormalised_Old",days[d],".txt"),header = T)
    tableExpressionOldCountNormalised<-merge(tableExpressionOldCountNormalised,expresstempo,by.x="gene",by.y="gene")
    expresstempo<-read.table(file = paste0(inDataNatmi,"CountNormalised_Young",days[d],".txt"),header = T)
    tableExpressionYoungCountNormalised<-merge(tableExpressionYoungCountNormalised,expresstempo,by.x="gene",by.y="gene")
  }
  tableExpressionYoungCountNormalised<- tableExpressionYoungCountNormalised %>% mutate(symbol=genes_df[match(tableExpressionYoungCountNormalised$gene,genes_df$Geneid),]$symbol)
  tableExpressionOldCountNormalised<- tableExpressionOldCountNormalised %>% mutate(symbol=genes_df[match(tableExpressionOldCountNormalised$gene,genes_df$Geneid),]$symbol)
  
  write_csv(tableExpressionYoungCountNormalised,paste0(inDataNatmi,concatExprAllReplicateYoung_CN))
  write_csv(tableExpressionOldCountNormalised,paste0(inDataNatmi,concatExprAllReplicateOld_CN))
  
  
  tableExpressionYoungMeanTPM<-tableExpressionYoungTPM %>% select(gene,symbol)
  tableExpressionOldMeanTPM<-tableExpressionOldTPM %>% select(gene,symbol)
  tableExpressionYoungMeanCN<-tableExpressionYoungCountNormalised %>% select(gene,symbol)
  tableExpressionOldMeanCN<-tableExpressionOldCountNormalised %>% select(gene,symbol)
  for (t in Typecellv){
    for (d in daysv[[t]]){
      tableExpressionYoungMeanTPM[[paste0("Young.Mean.TPM.",t,".",d)]]<- tableExpressionYoungTPM %>% select(contains(paste0("Young.",t,".",d))) %>% rowMeans() 
      tableExpressionOldMeanTPM[[paste0("Old.Mean.TPM.",t,".",d)]]<- tableExpressionOldTPM %>% select(contains(paste0("Old.",t,".",d))) %>% rowMeans() 
      tableExpressionYoungMeanCN[[paste0("Young.Mean.CN.",t,".",d)]]<- tableExpressionYoungCountNormalised %>% select(contains(paste0("Young.",t,".",d))) %>% rowMeans() 
      tableExpressionOldMeanCN[[paste0("Old.Mean.CN.",t,".",d)]]<- tableExpressionOldCountNormalised %>% select(contains(paste0("Old.",t,".",d))) %>% rowMeans() 
    }
  }
  for (d in days){
    tableExpressionYoungMeanTPM[[paste0("Young.Sum.TPM.",d)]]<- tableExpressionYoungMeanTPM %>% select(contains(d)) %>% rowSums() 
    tableExpressionOldMeanTPM[[paste0("Old.Sum.TPM.",d)]]<- tableExpressionOldMeanTPM %>% select(contains(d)) %>% rowSums() 
    tableExpressionYoungMeanCN[[paste0("Young.Sum.CN.",d)]]<- tableExpressionYoungMeanCN %>% select(contains(d)) %>% rowSums() 
    tableExpressionOldMeanCN[[paste0("Old.Sum.CN.",d)]]<- tableExpressionOldMeanCN %>% select(contains(d)) %>% rowSums() 
  }
  write_csv(tableExpressionYoungMeanTPM,paste0(inDataNatmi,ExprMeanTPMYoung_SumDay))
  write_csv(tableExpressionYoungMeanCN,paste0(inDataNatmi,ExpreMeanCNYoung_SumDay))
  write_csv(tableExpressionOldMeanTPM,paste0(inDataNatmi,ExprMeanTPMOld_SumDay))
  write_csv(tableExpressionOldMeanCN,paste0(inDataNatmi,ExpreMeanCNOld_SumDay))
  
} else{
  tableExpressionYoungMeanTPM<- read.csv(paste0(inDataNatmi,ExprMeanTPMYoung_SumDay))
  tableExpressionYoungMeanCN<-read.csv(paste0(inDataNatmi,ExpreMeanCNYoung_SumDay))
  tableExpressionOldMeanTPM<-read.csv(paste0(inDataNatmi,ExprMeanTPMOld_SumDay))
  tableExpressionOldMeanCN<-read.csv(paste0(inDataNatmi,ExpreMeanCNOld_SumDay))
}
  

#Concat NatmiOutFile
tableOld_CN=data.frame()
tableYoung_CN=data.frame()
tableDiff_CN=data.frame()
for (d in days ){
  tempoOld<-read.csv(paste0(odir,natmiOut_CN,"Old",d,"/",EdgesNatmi,"Edges.csv"))
  tempoOld<-tempoOld %>% mutate(day=rep(d,length(tempoOld$Ligand.symbol)))
  tableOld_CN<-rbind(tableOld_CN,tempoOld)
  tempoYoung<-read.csv(paste0(odir,natmiOut_CN,"Young",d,"/",EdgesNatmi,"Edges.csv"))
  tempoYoung<-tempoYoung %>% mutate(day=rep(d,length(tempoYoung$Ligand.symbol)))
  tableYoung_CN<-rbind(tableYoung_CN,tempoYoung)
  tempoDiff<-read.csv(paste0(odir,natmiOut_CN,"Diff_",d,"/",DeltaNatmi,"All_edges_mean.csv"))
  tempoDiff<-tempoDiff %>% mutate(day=rep(d,length(tempoDiff$Ligand.symbol)))
  tableDiff_CN<-rbind(tableDiff_CN,tempoDiff)
}
colnames(tableDiff_CN)<-str_replace_all(colnames(tableDiff_CN),".in.condition.1",".Young")
colnames(tableDiff_CN)<-str_replace_all(colnames(tableDiff_CN),".in.condition.2",".Old")

saveRDS(tableDiff_CN,paste0(odir,natmiOut_CN,"rds/TableDiff.rds"))
saveRDS(tableOld_CN,paste0(odir,natmiOut_CN,"rds/TableOld.rds"))
saveRDS(tableYoung_CN,paste0(odir,natmiOut_CN,"rds/TableYoung.rds"))

write.csv(tableDiff_CN,paste0(odir,natmiOut_CN,"csv/TableDiff.csv"))
write.csv(tableOld_CN,paste0(odir,natmiOut_CN,"csv/TableOld.csv"))
write.csv(tableYoung_CN,paste0(odir,natmiOut_CN,"csv/TableYoung.csv"))

tableOld_TPM=data.frame()
tableYoung_TPM=data.frame()
tableDiff_TPM=data.frame()
for (d in days ){
  tempoOld<-read.csv(paste0(odir,natmiOut_TPM,"Old",d,"/",EdgesNatmi,"Edges.csv"))
  tempoOld<-tempoOld %>% mutate(day=rep(d,length(tempoOld$Ligand.symbol)))
  tableOld_TPM<-rbind(tableOld_TPM,tempoOld)
  tempoYoung<-read.csv(paste0(odir,natmiOut_TPM,"Young",d,"/",EdgesNatmi,"Edges.csv"))
  tempoYoung<-tempoYoung %>% mutate(day=rep(d,length(tempoYoung$Ligand.symbol)))
  tableYoung_TPM<-rbind(tableYoung_TPM,tempoYoung)
  tempoDiff<-read.csv(paste0(odir,natmiOut_TPM,"Diff_",d,"/",DeltaNatmi,"All_edges_mean.csv"))
  tempoDiff<-tempoDiff %>% mutate(day=rep(d,length(tempoDiff$Ligand.symbol)))
  tableDiff_TPM<-rbind(tableDiff_TPM,tempoDiff)
}
colnames(tableDiff_TPM)<-str_replace_all(colnames(tableDiff_TPM),".in.condition.1",".Young")
colnames(tableDiff_TPM)<-str_replace_all(colnames(tableDiff_TPM),".in.condition.2",".Old")

saveRDS(tableDiff_TPM,paste0(odir,natmiOut_TPM,"rds/TableDiff.rds"))
saveRDS(tableOld_TPM,paste0(odir,natmiOut_TPM,"rds/TableOld.rds"))
saveRDS(tableYoung_TPM,paste0(odir,natmiOut_TPM,"rds/TableYoung.rds"))

write.csv(tableDiff_TPM,paste0(odir,natmiOut_TPM,"csv/TableDiff.csv"))
write.csv(tableOld_TPM,paste0(odir,natmiOut_TPM,"csv/TableOld.csv"))
write.csv(tableYoung_TPM,paste0(odir,natmiOut_TPM,"csv/TableYoung.csv"))


########################################################################################
# weight edge with ratio
# of expression Ligand on the sum of all ligand in couple with the receptor
# and ratio expression of Receptor on the sum af all receptor in couple with ligand
########################################################################################



Ligand_Receptors<-lapply(unique(tableDiff_TPM$Ligand.symbol), function(L) tableDiff_TPM %>% filter(Ligand.symbol==L) %>% select(Receptor.symbol) %>% unique() %>% unlist() %>% as.character())
names(Ligand_Receptors)<-unique(tableDiff_TPM$Ligand.symbol)
SumTPMReceptorY<-lapply(days, function(D) lapply(Ligand_Receptors, function(R) tableExpressionYoungMeanTPM %>% filter(symbol %in% R) %>% select(paste0("Young.Sum.TPM.",D)) %>% colSums() %>% unlist() %>% as.numeric() %>% round(digits = 8)))
names(SumTPMReceptorY)[]<-days
names(SumTPMReceptorY)[]<-days
SumTPMReceptorO<-lapply(days, function(D) lapply(Ligand_Receptors, function(R) tableExpressionOldMeanTPM %>% filter(symbol %in% R) %>% select(paste0("Old.Sum.TPM.",D)) %>% colSums() %>% unlist() %>% as.numeric() %>% round(digits = 8) ))
names(SumTPMReceptorO)[]<-days
names(SumTPMReceptorO)[]<-days

Receptor_Ligands<-lapply(unique(tableDiff_TPM$Receptor.symbol), function(R) tableDiff_TPM %>% filter(Receptor.symbol==R) %>% select(Ligand.symbol) %>% unique() %>% unlist() %>% as.character())
names(Receptor_Ligands)<-unique(tableDiff_TPM$Receptor.symbol)
SumTPMLigandY<-lapply(days, function(D) lapply(Receptor_Ligands, function(R) tableExpressionYoungMeanTPM %>% filter(symbol %in%  R) %>% select(paste0("Young.Sum.TPM.",D)) %>% colSums() %>% unlist() %>% as.numeric() %>% round(digits = 8)))
names(SumTPMLigandY)[]<-days
names(SumTPMLigandY)[]<-days
SumTPMLigandO<-lapply(days, function(D) lapply(Receptor_Ligands, function(R) tableExpressionOldMeanTPM %>% filter(symbol %in% R) %>% select(paste0("Old.Sum.TPM.",D)) %>% colSums() %>% unlist() %>% as.numeric() %>% round(digits = 8)))
names(SumTPMLigandO)[]<-days
names(SumTPMLigandO)[]<-days

tableCalcRatio<-tableDiff_TPM %>% select(Sending.cluster:Target.cluster,Ligand.expression.Young,Receptor.expression.Young,Ligand.expression.Old,Receptor.expression.Old,day) %>% filter(day=="D0")
tableCalcRatio <- tableCalcRatio %>% mutate()
tableCalcRatio<-tableCalcRatio %>%
  mutate( Ligands_SumYoungTPM=SumTPMLigandY[["D0"]][match(tableCalcRatio$Receptor.symbol,names(SumTPMLigandY[["D0"]]))]) %>%
  mutate(RatioLigandY= ifelse(is.infinite(as.numeric(Ligand.expression.Young)/as.numeric(Ligands_SumYoungTPM)),0.00001 ,as.numeric(Ligand.expression.Young)/as.numeric(Ligands_SumYoungTPM)))  %>%
  mutate(Receptors_SumYoungTPM=SumTPMReceptorY[["D0"]][match(tableCalcRatio$Ligand.symbol,names(SumTPMReceptorY[["D0"]]))]) %>%
  mutate(RatioReceptorY=ifelse(is.infinite(as.numeric(Receptor.expression.Young)/as.numeric(Receptors_SumYoungTPM)),0.00001,as.numeric(Receptor.expression.Young)/as.numeric(Receptors_SumYoungTPM))) %>%
  mutate(weightLigandYoung=RatioReceptorY*Ligand.expression.Young) %>%
  mutate(weightReceptorYoung=RatioLigandY*Receptor.expression.Young) %>%
  mutate(weightEdgeYoung=RatioLigandY*Receptor.expression.Young*RatioReceptorY*Ligand.expression.Young) %>%
  mutate( Ligands_SumOldTPM=SumTPMLigandO[["D0"]][match(tableCalcRatio$Receptor.symbol,names(SumTPMLigandO[["D0"]]))]) %>%
  mutate(RatioLigandO=ifelse(is.infinite(as.numeric(Ligand.expression.Old)/as.numeric(Ligands_SumOldTPM)),0.00001,as.numeric(Ligand.expression.Old)/as.numeric(Ligands_SumOldTPM))) %>%
  mutate( Receptors_SumOldTPM=SumTPMReceptorO[["D0"]][match(tableCalcRatio$Ligand.symbol,names(SumTPMReceptorO[["D0"]]))]) %>%
  mutate(RatioReceptorO=ifelse(is.infinite(as.numeric(Receptor.expression.Old)/as.numeric(Receptors_SumOldTPM)),0.00001,as.numeric(Receptor.expression.Old)/as.numeric(Receptors_SumOldTPM))) %>%
  mutate(weightLigandOld=RatioReceptorO*Ligand.expression.Old) %>%
  mutate(weightReceptorOld=RatioLigandO*Receptor.expression.Old) %>%
  mutate(weightEdgeOld=RatioLigandO*Receptor.expression.Old*RatioReceptorO*Ligand.expression.Old) %>%
  mutate(weightEdgeDiff=weightEdgeOld-weightEdgeYoung) %>%
  mutate(log2weightEdgeDiff=weightEdgeDiff/abs(weightEdgeDiff)*log2(abs(weightEdgeDiff)+1+10e-14))
for ( D in c("D2","D4","D7")){
  tableCalcRatiotempo<-tableDiff_TPM %>% select(Sending.cluster:Target.cluster,Ligand.expression.Young,Receptor.expression.Young,Ligand.expression.Old,Receptor.expression.Old,day) %>% filter(day==D)
  tableCalcRatiotempo<-tableCalcRatiotempo %>% mutate( Ligands_SumYoungTPM=SumTPMLigandY[[D]][match(tableCalcRatiotempo$Receptor.symbol,names(SumTPMLigandY[[D]]))]) %>%
    mutate(RatioLigandY=ifelse(is.infinite(as.numeric(Ligand.expression.Young)/as.numeric(Ligands_SumYoungTPM)),0.00001,as.numeric(Ligand.expression.Young)/as.numeric(Ligands_SumYoungTPM))) %>%
    mutate(Receptors_SumYoungTPM=SumTPMReceptorY[[D]][match(tableCalcRatiotempo$Ligand.symbol,names(SumTPMReceptorY[[D]]))]) %>%
    mutate(RatioReceptorY=ifelse(is.infinite(as.numeric(Receptor.expression.Young)/as.numeric(Receptors_SumYoungTPM)),0.00001,as.numeric(Receptor.expression.Young)/as.numeric(Receptors_SumYoungTPM))) %>%
    mutate(weightLigandYoung=RatioReceptorY*Ligand.expression.Young) %>%
    mutate(weightReceptorYoung=RatioLigandY*Receptor.expression.Young) %>%
    mutate(weightEdgeYoung=RatioLigandY*Receptor.expression.Young*RatioReceptorY*Ligand.expression.Young) %>%
    mutate( Ligands_SumOldTPM=SumTPMLigandO[[D]][match(tableCalcRatiotempo$Receptor.symbol,names(SumTPMLigandO[[D]]))]) %>%
    mutate(RatioLigandO=ifelse(is.infinite(as.numeric(Ligand.expression.Old)/as.numeric(Ligands_SumOldTPM)),0.00001,as.numeric(Ligand.expression.Old)/as.numeric(Ligands_SumOldTPM))) %>%
    mutate( Receptors_SumOldTPM=SumTPMReceptorO[[D]][match(tableCalcRatiotempo$Ligand.symbol,names(SumTPMReceptorO[[D]]))]) %>%
    mutate(RatioReceptorO=ifelse(is.infinite(as.numeric(Receptor.expression.Old)/as.numeric(Receptors_SumOldTPM)),0.00001,as.numeric(Receptor.expression.Old)/as.numeric(Receptors_SumOldTPM))) %>%
    mutate(weightLigandOld=RatioReceptorO*Ligand.expression.Old) %>%
    mutate(weightReceptorOld=RatioLigandO*Receptor.expression.Old) %>%
    mutate(weightEdgeOld=RatioLigandO*Receptor.expression.Old*RatioReceptorO*Ligand.expression.Old) %>%
    mutate(weightEdgeDiff=weightEdgeOld-weightEdgeYoung)  %>%
    mutate(log2weightEdgeDiff=weightEdgeDiff/abs(weightEdgeDiff)*log2(abs(weightEdgeDiff)+1+10e-14))
  tableCalcRatio<- rbind(tableCalcRatio,tableCalcRatiotempo)
}

#Explore Weight Receptor, Ligand, Edge  Young Old and Diff to find threshold Weight Edge Diff OY -> extract a list of RL 
## Explore Weight Receptor, Ligand, Edge  Young Old and Diff

tableCalcRatio2<- tableCalcRatio %>% pivot_longer(cols = starts_with("weight"), names_to = "LRED" ,  names_prefix = "weight", "values_to"= "weight")

ggplot(tableCalcRatio2 , aes(weight/abs(weight)*log2(abs(weight))+1+10e-14) ) + geom_histogram()+ facet_wrap(~LRED) 

tableCalcRatio3<- tableCalcRatio %>% pivot_longer(cols = starts_with("Ratio"), names_to = "RLYO" ,  names_prefix = "ratio", "values_to"= "Ratio")

ggplot(tableCalcRatio3 , aes(Ratio)) + geom_histogram()+ facet_wrap(~RLYO) 


#preliminary exploration of the list with upset plot of curated list on days



png(paste0(odir,natmiOut_TPM,"analysis/log2weightEdgeDiff_distribution.png"))
layout(mat = matrix(c(1,2),2,1, byrow=TRUE),  height = c(1,8))
# Draw the boxplot and the histogram 
par(mar=c(0, 3.1, 1.1, 2.1))
boxplot(tableCalcRatio$log2weightEdgeDiff , horizontal=TRUE , xaxt="n" , col=rgb(0.8,0.8,0,0.5) , frame=F)
par(mar=c(4, 3.1, 1.1, 2.1))
hist(tableCalcRatio$log2weightEdgeDiff , breaks=40 , col="#69b3a2" , border=F , main="" , xlab="log2weightEdgeDiff") 
dev.off()
thresholdDiffPos=10
thresholdDiffNeg=-10

tableCalcRatioTopUP <- tableCalcRatio %>% filter(weightEdgeDiff>0) %>% mutate(RLday=paste0(Receptor.symbol,"_",Ligand.symbol,'_',day)) %>% arrange(desc(log2weightEdgeDiff)) %>% head(n=20L)
tableCalcRatioTopDown <- tableCalcRatio %>% filter(weightEdgeDiff<0)  %>% mutate(RLday=paste0(Receptor.symbol,"_",Ligand.symbol,'_',day)) %>% arrange(log2weightEdgeDiff) %>% head(n=20L)

UmapWeightEdge <- ggplot(tableCalcRatio, aes(x=log2(weightEdgeYoung+1+10e-14), y=log2(weightEdgeOld+1+10e-14), color=log2weightEdgeDiff),alpha=0.5) +
  geom_point() +
  geom_text_repel(data= rbind(tableCalcRatioTopUP,tableCalcRatioTopDown),
                  aes(label=RLday),col="black", size=3, segment.size = .1, force=2, force_pull = 2, max.overlaps=40)+
  xlab("log2weightEdgeYoung")+
  ylab("log2weightEdgeOld")+
  theme_ipsum()+
  theme(legend.title = element_text(size=8))+
  ggtitle("Distribution of weight edge of receptor ligand young on old")+
  scale_color_gradient2(name="WeightEdgeDiffOY",midpoint=0,  low="blue", mid="white",
                        high="red")

UmapWeightEdge_Hist <- ggMarginal(UmapWeightEdge, type="histogram")
UmapWeightEdge_Dens <- ggMarginal(UmapWeightEdge, type="density")
UmapWeightEdge_Box <- ggMarginal(UmapWeightEdge, type="boxplot")

tableCalcRatio2<-tableCalcRatio %>% filter(log2weightEdgeDiff > thresholdDiffPos | log2weightEdgeDiff < thresholdDiffNeg)
UmapWeightEdge2 <- ggplot(tableCalcRatio2, aes(x=log2(weightEdgeYoung+1+10e-14), y=log2(weightEdgeOld+1+10e-14), color=log2weightEdgeDiff),alpha=0.5) +
  geom_point() +
  geom_text_repel(data= rbind(tableCalcRatioTopUP,tableCalcRatioTopDown),
                  aes(label=RLday),col="black", size=3, segment.size = .1, force=2, force_pull = 2, max.overlaps=40)+
  xlab("log2weightEdgeYoung")+
  ylab("log2weightEdgeOld")+
  theme_ipsum()+
  theme(legend.title = element_text(size=8))+
  ggtitle("Distribution of weight edge of receptor ligand young on old")+
  scale_color_gradient2(name="WeightEdgeDiffOY",midpoint=0,  low="blue", mid="grey88",
                        high="red")

UmapWeightEdge_Hist2 <- ggMarginal(UmapWeightEdge2, type="histogram")
UmapWeightEdge_Dens2 <- ggMarginal(UmapWeightEdge2, type="density")
UmapWeightEdge_Box2 <- ggMarginal(UmapWeightEdge2, type="boxplot")


tableDiff_TPMTopUp <- tableDiff_TPM %>% filter(Delta.edge.expression.weight>0) %>% mutate(Log2weightEdgeDiff=log2(abs(Delta.edge.expression.weight)+1+10e-14)) %>% mutate(RL=paste0(Receptor.symbol,"_",Ligand.symbol)) %>% arrange(desc(Log2weightEdgeDiff)) %>% head(n=25L)
tableDiff_TPMTopDown <- tableDiff_TPM %>% filter(Delta.edge.expression.weight<0) %>% mutate(Log2weightEdgeDiff=log2(abs(Delta.edge.expression.weight)+1+10e-14)) %>% mutate(RL=paste0(Receptor.symbol,"_",Ligand.symbol)) %>% arrange(desc(Log2weightEdgeDiff)) %>% head(n=25L)

r <- ggplot(tableDiff_TPM, aes(x=log2(Edge.expression.weight.Young+1+10e-14), y=log2(Edge.expression.weight.Old+1+10e-14), color=Delta.edge.expression.weight/abs(Delta.edge.expression.weight)*log2(abs(Delta.edge.expression.weight))+1+10e-14),alpha(0.5)) +
  geom_point() +
  scale_color_gradient2(midpoint=0,  low="blue", mid="white",
                        high="red")+
  theme(legend.position="none") +
  geom_text_repel(data= rbind(tableDiff_TPMTopUp,tableDiff_TPMTopDown),
                  aes(label=RL),col="black", size=3, segment.size = .1, force=2, force_pull = 2, max.overlaps=15)

r1 <- ggMarginal(p, type="histogram")

write_csv(tableCalcRatio,paste0(odir,natmiOut_TPM,"analysis/tableDiffWeightEdge.csv"))
save_plot(paste0(odir,natmiOut_TPM,"analysis/WeightEdgeYoungOnWeightEdgeOld.png"),p1,base_height = 8, base_width = 15)

### Extract list of Receptor ligand DE , and start analysis with raw count normalised
## Integrate DEG info, add the ligand majoritaire et the receptor majoritaire of couple
## 


######################
# Add information DEG
######################

RLsigni<-data.frame()
RLstat<-data.frame()


TypeDay<-lapply( names(daysv) , function(x) paste0(x,".",daysv[[x]]) ) %>% unlist()
tableExpressionYoungMeanReplicate<-tableExpressionYoungCountNormalised %>% select(gene,symbol)
tableExpressionYoungMeanReplicate= cbind(tableExpressionYoungMeanReplicate,lapply(TypeDay, function(x) tableExpressionYoungCountNormalised %>% select(contains(x)) %>% rowMeans()))
colnames(tableExpressionYoungMeanReplicate)<-c("gene","symbol",paste0("Young.",TypeDay))
tableExpressionYoungMeanReplicate <- tableExpressionYoungMeanReplicate %>% mutate(Young.Sum.CountNormalised.D0=tableExpressionYoungMeanReplicate %>% select(contains(".D0")) %>% rowSums()) %>%
  mutate(Young.Sum.CountNormalised.D2=tableExpressionYoungMeanReplicate %>% select(contains(".D2")) %>% rowSums()) %>%
  mutate(Young.Sum.CountNormalised.D4=tableExpressionYoungMeanReplicate %>% select(contains(".D4")) %>% rowSums()) %>%
  mutate(Young.Sum.CountNormalised.D7=tableExpressionYoungMeanReplicate %>% select(contains(".D7")) %>% rowSums())
  
tableExpressionOldMeanReplicate<-tableExpressionOldCountNormalised %>% select(gene,symbol)
tableExpressionOldMeanReplicate= cbind(tableExpressionOldMeanReplicate,lapply(TypeDay, function(x) tableExpressionOldCountNormalised %>% select(contains(x)) %>% rowMeans()))
colnames(tableExpressionOldMeanReplicate)<-c("gene","symbol",paste0("Young.",TypeDay))
tableExpressionOldMeanReplicate <- tableExpressionOldMeanReplicate %>% mutate(Old.Sum.CountNormalised.D0=tableExpressionOldMeanReplicate %>% select(contains(".D0")) %>% rowSums()) %>%
  mutate(Old.Sum.CountNormalised.D2=tableExpressionOldMeanReplicate %>% select(contains(".D2")) %>% rowSums()) %>%
  mutate(Old.Sum.CountNormalised.D4=tableExpressionOldMeanReplicate %>% select(contains(".D4")) %>% rowSums()) %>%
  mutate(Old.Sum.CountNormalised.D7=tableExpressionOldMeanReplicate %>% select(contains(".D7")) %>% rowSums())


RLInfoSupp<- tableDiff %>%  mutate(DayTypeLigand=paste0(Ligand.symbol,"_",day,"_",Sending.cluster)) %>% mutate(DayTypeReceptor=paste0(Receptor.symbol,"_",day,"_",Target.cluster))
RLInfoSupp<- RLInfoSupp %>% mutate(Ligand.DEG.padj=fullDEsta[match(DayTypeLigand,fullDEsta$GeneDayType),]$padj) %>%
  mutate(Ligand.DEG.log2FoldChange=fullDEsta[match(DayTypeLigand,fullDEsta$GeneDayType),]$log2FoldChange) %>%
  mutate(Receptor.DEG.padj=fullDEsta[match(DayTypeReceptor,fullDEsta$GeneDayType),]$padj) %>%
  mutate(Receptor.DEG.log2FoldChange=fullDEsta[match(DayTypeReceptor,fullDEsta$GeneDayType),]$log2FoldChange) %>%
  mutate(Ligand.DEG=ifelse(as.numeric(Ligand.DEG.padj) < 0.05, "DEG", "noDEG")) %>%
  mutate(Receptor.DEG=ifelse(as.numeric(Receptor.DEG.padj) < 0.05, "DEG", "noDEG")) %>%
  mutate(RL=paste0(Receptor.symbol,"_",Ligand.symbol)) %>%
  subset(select=c(RL,day,Sending.cluster,Ligand.symbol,Ligand.DEG,Ligand.DEG.padj,Ligand.DEG.log2FoldChange,Target.cluster,Receptor.symbol,Receptor.DEG,Receptor.DEG.padj,Receptor.DEG.log2FoldChange,Ligand.expressing.cells.Young:Log2.transformed.fold.change.of.edge.specificity.weight)) 


RLDiffSumExpress<-data.frame( RL= rep(unique(RLInfoSupp$RL),4),Day=c(rep("D0",length(unique(RLInfoSupp$RL))),rep("D2",length(unique(RLInfoSupp$RL))),rep("D4",length(unique(RLInfoSupp$RL))),rep("D7",length(unique(RLInfoSupp$RL))))) %>%
                 mutate(Receptor.symbol= str_split_fixed(RL,"_",2)[,1] , Ligand.symbol =str_split_fixed(RL,"_",2)[,2])

getNormalisedExpress<-function(r,tableExpress,age,col){
  d=RLDiffSumExpress[r,]$Day 
  symbolReceptor=RLDiffSumExpress[r,col] 
  expr<-tableExpress %>% filter(symbol==symbolReceptor) %>% select(paste0(age,".Sum.CountNormalised.",d))  %>% unlist() %>% as.numeric()
  return(expr)}
RLDiffSumExpress<- RLDiffSumExpress %>% mutate(Receptor.Young.Count.Normalised= tableExpressionYoungMeanReplicate[match(RLDiffSumExpress$Receptor.symbol,tableExpressionYoungMeanReplicate$symbol),] %>% select(contains("Young.Sum.N")) lapply(rownames(RLDiffSumExpress), function(r) getNormalisedExpress(r,tableExpressionYoungMeanReplicate,"Young",3) )) %>%
  mutate(Receptor.Old.Count.Normalised=lapply(rownames(RLDiffSumExpress), function(r) getNormalisedExpress(r,tableExpressionOldMeanReplicate,"Old",3) )) %>%
  mutate(Ligand.Young.Count.Normalised=lapply(rownames(RLDiffSumExpress), function(r) getNormalisedExpress(r,tableExpressionYoungMeanReplicate,"Young",4) )) %>%
  mutate(Ligand.Old.Count.Normalised=lapply(rownames(RLDiffSumExpress), function(r) getNormalisedExpress(r,tableExpressionOldMeanReplicate,"Old",4) ))

RLDiffSumExpress<- RLDiffSumExpress %>% mutate(Nb.Cell.Express.Ligand.Young=  )



###################################################################
# BMP pathways LR analysis
###################################################################
BMPanalysisLR<-RLInfoSupp %>% filter(str_detect(Ligand.symbol,"Bmp")|
  Receptor.symbol==str_detect(Receptor.symbol,"Bmp"))

BMPanalysisLRsigni<-RLsigni %>% filter(str_detect(Ligand.symbol,"Bmp")|Receptor.symbol==str_detect(Receptor.symbol,"Bmp"))%>% 
  mutate(SumDetectionRate=Ligand.detection.rate.Young+Receptor.detection.rate.Young+Ligand.detection.rate.Old+Receptor.detection.rate.Old) %>%
  filter(SumDetectionRate>=2)


write.csv(BMPanalysisLR,paste0(odir,"BMP_receptor_ligand.csv"))
write.csv(BMPanalysisLRsigni,paste0(odir,"BMP_receptor_ligand_DEG.csv"))
write.csv(DEGinfo %>% filter(str_detect(symbol,"Smad")),paste0(odir,"Smad_DEG.csv"))
BMPanalysisLRsigni<-BMPanalysisLRsigni %>% mutate(SumLRFoldchange=abs(as.numeric(Ligand.DEG.log2FoldChange))+abs(as.numeric(Receptor.DEG.log2FoldChange)))



###################################################################
# plot network
###################################################################
read.graph(file=past
           ("/home/bioinfo/BulkAnalysis_plusNetwork/AgeingMice_shinyApp/graphobjs_copy", age,"_",currentday$x,"_igraph_unfi.ml"), 
           format="graphml")
for (age in c("Young","Old")){
  V_vec[[age]] =
    igraph::as_data_frame(igElems_list$x[[age]], 
                          "vertices")[["_nx_name"]]}
doinduced <- function(g, interestvec, orderinp){
  vertex_attr(g)$numid <- vertex_attr(g)$id
  vertex_attr(g)$uniname <- vertex_attr(g)[["_nx_name"]]
  selnodes = V(g)[uniname %in% interestvec]
  selegoV <- ego(g, order=orderinp, nodes = selnodes, mode = "all", mindist = 0)
  # turn the returned list of igraph.vs objects into a graph
  selegoG <- induced_subgraph(g,unlist(selegoV))
  return(selegoG)
}
doinduced(
  igElems_list$x[["Young"]], pickneigh[["Young"]], orderinp = input$NEIGH )
print(head(igraph::as_data_frame(here.subgr[["Young"]], "vertices"),1))
subnet.y$nodes <- massageDATA(here.subgr[["Young"]])$nodes
subnet.y$edges <- massageDATA(here.subgr[["Young"]])$edges
print("ok subnet young")
output$labtextyoung <- renderText({paste("YOUNG day ", currentday$x)})
output$young <- renderVisNetwork({
  req(subnet.y$edges)
  print(unique(subnet.y$nodes$groupname))
  
  netout <- visNetwork(subnet.y$nodes,
                       subnet.y$edges,
                       physics = FALSE) %>%
    visEdges(shadow=T, arrows="to", physics=FALSE) %>%
    visInteraction(navigationButtons = TRUE)
  netout
}) 
massageDATA <- function(myg){
  data <- toVisNetworkData(myg)
  data$nodes$label = data$nodes$genesym
  data$nodes$value = data$nodes$specificity * 10  # !! $averagexp changed!
  data$nodes$groupname = data$nodes$celltype
  data$edges$width = data$edges$weight * 10
  return(data)
}