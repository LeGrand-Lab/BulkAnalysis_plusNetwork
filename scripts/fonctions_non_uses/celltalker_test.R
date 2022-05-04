

suppressMessages({
  library(Seurat)
  library(celltalker)
})



library(tidyverse)
library(dplyr)
library(ggplot2)
library(stringr)

setwd("~/BulkAnalysis_plusNetwork/")
prefil_cou <- "data/prefiltered_counts.rds"
metadata.rds <- "data/metadata.rds"
genes_df <- read.table("data/genesinfo.csv", sep="\t", header=T)
dataGene = "exam_INTER_conditions/static/rds/"
DEGfile =  "TableDynamicUpDownDEG.rds"

fmat <- readRDS(prefil_cou)
rownames(fmat)<-genes_df[match(rownames(fmat), genes_df$Geneid),]$symbol
extractGeneDuplicate <- rownames(fmat[duplicated(rownames(fmat)),])
design<-readRDS(metadata.rds)
design <- design %>% mutate(type2=str_replace(type, "M1|M2|Neutro","cellImmun" ))
fmat<-fmat[!rownames(fmat) %in% extractGeneDuplicate,]
head(fmat)
colnames(fmat)<-str_replace_all(colnames(fmat), "_" , "-")

# Load DataBase receptor Ligand
load(file="/home/bioinfo/prog_bio/celltalker-master/data/ramilowski_pairs.rda")
CONNECTOMEDB<- read.table("/home/bioinfo/NATMI/lrdbs/lrc2a_v2.txt", header = T)
# Convert human gene symbol in mus musculus gene symbol
homology<-read.table("/home/bioinfo/NATMI/homology/homologene_9606_10090.data")
GENEID_9606<-homology %>% filter(V2=="9606") %>% select(V1,V4) 
GENEID_10090<-homology %>% filter(V2=="10090") %>% select(V1,V4) 

correspondance_9606_10090 <- inner_join(GENEID_9606,GENEID_10090, by="V1") 
colnames(correspondance_9606_10090) <- c("ID","9606","10090") 
correspondance_9606_10090 <- correspondance_9606_10090 %>% select("9606","10090")
write_csv(correspondance_9606_10090,"/home/bioinfo/NATMI/homology/correspondance_9606_10090.csv" )

CONNECTOMEDB_tempo<-merge(CONNECTOMEDB,correspondance_9606_10090 ,by.x="ligand", by.y="9606")
CONNECTOMEDB_tempo<-merge(CONNECTOMEDB_tempo,correspondance_9606_10090 ,by.x="receptor", by.y="9606")

write_csv(CONNECTOMEDB_tempo,"/home/bioinfo/NATMI/homology/CONNECTOMEDN_correspondance_9606_10090.csv" )
CONNECTOMEDB_mmu<-data.frame("ligand"=CONNECTOMEDB_tempo$`10090.x`,"receptor"=CONNECTOMEDB_tempo$`10090.y`,"pair"=paste0(CONNECTOMEDB_tempo$`10090.x`,"_",CONNECTOMEDB_tempo$`10090.y`) )

for (d in unique(design$time ) ){
  fmat <- as.data.frame(fmat) %>% select(contains(d))
  design_d <- design %>% filter(time %in% d)
  metadata <- data.frame("sample.id"=paste0(design_d$age,".",sapply(strsplit(design_d$newname,split="_"),function(x) x[2])),"sample.type"=design_d$age,"cell.type"=design_d$type2)
  rownames(metadata)<-colnames(fmat)
  
  ser.obj <- CreateSeuratObject(counts=fmat,meta.data=metadata)
  
  #Standard Seurat workflow
  
  ser.obj <- NormalizeData(ser.obj )
  ser.obj <- FindVariableFeatures(ser.obj)
  ser.obj <- ScaleData(ser.obj)
  
  #ser.obj <- RunPCA(ser.obj)
  #ElbowPlot(ser.obj)
  #ser.obj <- RunUMAP(ser.obj,reduction="pca",dims=1:15)
  #ser.obj <- FindNeighbors(ser.obj,reduction="pca",dims=1:15)
  #ser.obj <- FindClusters(ser.obj,resolution=0.5)
  #p1 <- DimPlot(ser.obj,reduction="umap",group.by="age") 
  #p2 <- DimPlot(ser.obj,reduction="umap",group.by="day" )
  #p3 <- DimPlot(ser.obj,reduction="umap",group.by="cell.type")
  #p4 <- DimPlot(ser.obj,reduction="umap",group.by="RNA_snn_res.0.5",label=T) + NoLegend()
  #cowplot::plot_grid(p1,p2,p3,p4)
  
  ligsCDB <- as.character(unique(CONNECTOMEDB_mmu$ligand))
  recsCDB <- as.character(unique(CONNECTOMEDB_mmu$receptor))
  
  ligs.present <- rownames(ser.obj)[rownames(ser.obj) %in% ligsCDB]
  recs.present <- rownames(ser.obj)[rownames(ser.obj) %in% recsCDB]
  
  genes.to.use <- union(ligs.present,recs.present)
  
  
  Idents(ser.obj) <- "cell.type"
  
  DEGfile <- read.csv("/home/bioinfo/BulkAnalysis_plusNetwork/exam_INTER_conditions/static/csv/TableDynamicDEG.csv")
  ligs.recs.use <- DEGfile %>% filter(day==d, symbol %in% genes.to.use) %>% select(symbol) %>% unlist() %>% as.character() %>% unique()
  length(ligs.recs.use)
  
  interactions.forward1 <- CONNECTOMEDB_mmu[as.character(CONNECTOMEDB_mmu$ligand) %in% genes.to.use,]
  interactions.forward2 <- CONNECTOMEDB_mmu[as.character(CONNECTOMEDB_mmu$receptor) %in% genes.to.use,]
  interact.for <- rbind(interactions.forward1,interactions.forward2)
  dim(interact.for)             
  
  expr.mat <- GetAssayData(ser.obj,slot="counts")
  defined.clusters <- ser.obj@meta.data$cell.type
  defined.groups <- ser.obj@meta.data$sample.type
  defined.replicates <- ser.obj@meta.data$sample.id
  
  reshaped.matrices <- celltalker::reshape_matrices(count.matrix=expr.mat,clusters=as.factor(defined.clusters),groups=as.factor(defined.groups),replicates=as.factor(defined.replicates),ligands.and.receptors=interact.for)
  
  #Check out the hierarchy of the tibble
  reshaped.matrices
  unnest(reshaped.matrices,cols="samples")
  names(pull(unnest(reshaped.matrices,cols="samples"))[[1]])
  consistent.lig.recs <- create_lig_rec_tib(exp.tib=reshaped.matrices,clusters=defined.clusters,groups=defined.groups,replicates=defined.replicates,cells.reqd=2,freq.pos.reqd=0.3,ligands.and.receptors=interact.for)
  
  unnest(consistent.lig.recs[1,2],cols="lig.rec.exp")
  pull(unnest(consistent.lig.recs[1,2],cols="lig.rec.exp")[1,2])[[1]]
  
  put.int <- putative_interactions(ligand.receptor.tibble=consistent.lig.recs,clusters=defined.clusters,groups=defined.groups,freq.group.in.cluster=0.005,ligands.and.receptors=interact.for)
  unique.ints <- unique_interactions(put.int,group1="Young",group2="Old",interact.for) 
  
  pbmc.to.plot <- pull(unique.ints[1,2])[[1]]
  for.circos.pbmc <- pull(put.int[1,2])[[1]][pbmc.to.plot]
  
  circos_plot(interactions=for.circos.pbmc,clusters=defined.clusters)
}


