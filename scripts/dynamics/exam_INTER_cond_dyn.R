## 
# Detect Differential Expression between 2 times point consecutive for young, old and young vs old 
# output: all results into 'exam_INTER_conditions/dynamic/'
# --
# Joha GL - MoulleP
##
library(dplyr)
library(tidyverse)
library(cowplot)
library(forcats)
library(RColorBrewer)
display.brewer.all(colorblindFriendly = T)
library(pheatmap)
library(DESeq2)
library(ggsci) # publishing palettes
library(gridExtra)
library(reshape2)
library(ggthemes)
library(scales) # 'viridis_pal' et autres pal
library(ggrepel) # for labels to points
library(grid)
#remotes::install_github('coolbutuseless/facetious')
library(patchwork)
library(facetious)
library("apeglm") # BiocManager::install("apeglm")
library("BiocParallel")
library("ashr") # install.packages("ashr")
library(wesanderson)
register(MulticoreParam(4)) # TODO:  set n of cores depending of available
# devtools::install_github("vqv/ggbipplot")
# install.packages("glmpca")

##################################
#### Path files to load or created
##################################

#Directories
setwd("~/BulkAnalysis_plusNetwork/")
inData<- "data/"
NormData <- "data/CountNormalised/"
odir <- "exam_INTER_conditions/"
dynamicD =  paste0(odir,"dynamic/")
dynamicYD =paste0(dynamicD,"Young/")
dynamicOD =paste0(dynamicD,"Old/")

listDirectorie<-c(inData,NormData,odir,dynamicYD,dynamicOD,dynamicD)

for ( directories in listDirectorie){
  if (dir.exists(directories) == F) {
    dir.create(directories)
  }
}
listInterCondtionDirectories<-c(dynamicYD,dynamicOD,dynamicD)
for ( sd in c("TableDEG/","PlotsDEG/","analysis/") ){
  for (d in listInterCondtionDirectories ){
    if (dir.exists(paste0(d,sd)) == F) {
      dir.create(paste0(d,sd))
    }
  }
}

#Loads
genes_df <- read.table("data/genesinfo.csv", sep="\t", header=T)
prefil_cou <- paste0(inData,"prefiltered_counts.rds")
metadata.rds <- paste0(inData,"metadata.rds")

#Created
## static|dynamic,DEG_table,[Young_|Old_|][static|dynamic]_[full|solftfilter|hardfilter].rds")
## full -> all deseq2 analysis
## solt filter -> padj < 0.5 and abs( lfc ) > 0.1
## final filter -> padj < .05 and abs ( lfc ) >= 1.2
##Volvanoplot full table, text hard filter
# Name : [Young_|Old_|]volcano_static|dynamic

DEG_table = "TableDEG/DEG_table_"
volcanoneeded <- TRUE
orderTypecelldyn=c("ECs","FAPs","MuSCs","Inflammatory-Mac","Resolving-Mac")

colorsTypedyn=c("#10b387ff","#3d85c6ff","#b171f1ff","#ff9900ff","#cc0000ff")
names(colorsTypedyn)=orderTypecelldyn
colorsTimeCondition = c(brewer.pal(9,"RdPu")[3],brewer.pal(9,"RdPu")[5],brewer.pal(9,"RdPu")[7])


## Load files
fmat <- readRDS(prefil_cou)
metadata <- readRDS(metadata.rds)
metadata <- metadata %>% mutate(day=time) %>% mutate(daytype=paste0(day,".",type)) 
keep <- apply(fmat, 1, function(row) ifelse(length( which( row >= 5)) >= 3, TRUE, FALSE) )
fmat <- fmat[keep,]


#####
# DE analysis
# dynamic analyse for old and young
# by cell-type between t-1 and t
#####

##  with 'subset' s   on same original object:
#  ===========================================================================
timeline<-c("D0","D2","D4","D7")
typesv<-lapply( unique(metadata$day), function(d) unique(metadata[metadata$day == d,]$type) )
names(typesv)<-unique(metadata$day)
#All replicats
countNormalised2<-data.frame("gene"=rownames(fmat))
#Mean count t and t+1 on two line
countNormalisedtable<-data_frame(geneID = character(),age=character(),celltype=character(),TimeCondtion=character(),day=character(),meancountnormalised=character())
#Mean count t and t+1 on one line
countNormalisedtable2<-data_frame(geneID = character(),age=character(),celltype=character(),TimeCondtion=character(),meancountnormalised_t=character(),"meancountnormalised_t+1"=character())
for ( A in unique(metadata$age)){
  print(A)
  fmat_A<-as.data.frame(fmat) %>% dplyr::select( contains(A))
  metadata_A<-metadata  %>% dplyr::filter( age == A)
  sortie_full <- data.frame("baseMean"= double(), "log2FoldChange"= double(), 
                            "lfcSE"= double(),  "stat" = double(), 
                            "pvalue" = double(), "padj" = double(),
                            "id" = character(), "age" = character(),"TimeCondtion" = character(), "type"=character())
  out_softfilt<-data.frame()
  for ( t in 1:(length(timeline)-1)){
    print(t)
    D<-timeline[t]
    D2<-timeline[t+1]
    typeShared<-intersect(typesv[[D]], typesv[[D2]])
    fmat_AD<- fmat_A %>% dplyr::select(contains(c(D,D2))) %>% dplyr::select(contains(typeShared))  %>% as.matrix()
    metadata_AD <- metadata_A %>% dplyr::filter( day==D|day==D2, type %in% typeShared)    %>% mutate(TimeCondtion=ifelse(day == D , "TimeRef", "TimePlusOne")) %>% arrange(factor(newname,levels = colnames(fmat_AD)))
    
    deseqdataset.dynamic <- DESeqDataSetFromMatrix(countData = fmat_AD,
                                                   colData = metadata_AD,
                                                   design= ~ TimeCondtion)
    
    
    
    # set "Young" as reference
    colData(deseqdataset.dynamic)$TimeCondtion <- factor(colData(deseqdataset.dynamic)$TimeCondtion,
                                                         levels=c("TimeRef", "TimePlusOne"))
    
    
    for (t in typeShared){
      print(paste0(D,' vs ',D2,' in ', t))
      # do subsets as needed + deseq analysis
      deseqdataset.dynamic.daytype <- deseqdataset.dynamic[,deseqdataset.dynamic$type == t ]
      deseqdataset.dynamic.daytype <- DESeq(deseqdataset.dynamic.daytype, test="Wald", full=~TimeCondtion)
      
      countNormalised<-counts(deseqdataset.dynamic.daytype , normalize =T)
      countNormalised2<-cbind(countNormalised,countNormalised2)
      
      countNormalisedtable_tempo= data.frame(geneID=c(rownames(countNormalised),rownames(countNormalised)),
                                             age=c(rep(A, dim(countNormalised)[1]*2)),
                                             celltype=rep(t, dim(countNormalised)[1]*2),
                                             TimeCondtion=rep(paste0(D,'_',D2), dim(countNormalised)[1]*2),
                                             day=c(rep(D, dim(countNormalised)[1]),rep(D2, dim(countNormalised)[1])),
                                             meancountnormalised=c(rowMeans(countNormalised[,str_detect(colnames(countNormalised),D)]),rowMeans(countNormalised[,str_detect(colnames(countNormalised),D2)])))
      countNormalisedtable_tempo2= data.frame(geneID=rownames(countNormalised),
                                              age=rep(A, dim(countNormalised)[1]),
                                              celltype=rep(t, dim(countNormalised)[1]),
                                              TimeCondtion=rep(paste0(D,'_',D2), dim(countNormalised)[1]),
                                              meancountnormalised_t=rowMeans(countNormalised[,str_detect(colnames(countNormalised),D)]),
                                              "meancountnormalised_t+1"=rowMeans(countNormalised[,str_detect(colnames(countNormalised),D2)]))
      countNormalisedtable<- rbind(countNormalisedtable,countNormalisedtable_tempo)
      countNormalisedtable2<- rbind(countNormalisedtable2,countNormalisedtable_tempo2)
      # Add Shrink log2 fold changes
      deres <- lfcShrink(deseqdataset.dynamic.daytype, coef="TimeCondtion_TimePlusOne_vs_TimeRef", type="apeglm", parallel=T)
      deres$id <- rownames(deres)
      deres$TimeCondtion = paste0(D,"_",D2)
      deres$type = t
      deres$age= A
      
      sortie_full <- rbind(sortie_full, as_tibble(deres))
      
      deres_f <- as_tibble(deres) %>% dplyr::filter(padj <= 0.5 & abs(log2FoldChange) >= 0.1)
      out_softfilt <- rbind(out_softfilt,deres_f)
      print(head(sortie_full))
    }
  }
  
  # complete tables with gene_symbols & save table
  sortie_full$symbol <- genes_df[match(sortie_full$id, genes_df$Geneid),]$symbol
  sortie_full<- sortie_full %>% dplyr::filter(baseMean!=0) 
  saveRDS(sortie_full, paste0(dynamicD,A,"/", DEG_table, A,"_dynamic_full.rds")) 
  write.table(sortie_full, paste0(dynamicD,A,"/",DEG_table, A,"_dynamic_full.csv"), sep="\t", col.names=T, row.names = F) 
  out_softfilt <- sortie_full %>% dplyr::filter(padj <= 0.5 & abs(log2FoldChange) >= 0.1)
  write.table(out_softfilt, paste0(dynamicD,A,"/",DEG_table, A,"_dynamic_softfilter.csv"), sep="\t", col.names=T, row.names = F) 
  out_hardfilt = out_softfilt %>% dplyr::filter(padj <= 0.05 & abs(log2FoldChange) >= 1.2)
  write.table(out_hardfilt, paste0(dynamicD,A,"/",DEG_table, A,"_dynamic_hardfilter.csv"), sep="\t", col.names=T, row.names = F)
}
sortie_full_young<-readRDS( paste0(dynamicD,"Young/", DEG_table,"Young_dynamic_full.rds"))
sortie_full_old<-readRDS( paste0(dynamicD,"Old/", DEG_table,"Old_dynamic_full.rds"))
sortie_full_young<- sortie_full_young %>% mutate(ID=paste0(symbol,type,TimeCondtion,"Young"))
sortie_full_old<-sortie_full_old %>% mutate(ID=paste0(symbol,type,TimeCondtion,"Old"))
sortie_full<-rbind(sortie_full_old,sortie_full_young)


countNormalisedtable$symbol<- genes_df[match(countNormalisedtable$geneID, genes_df$Geneid),]$symbol
countNormalisedtable2$symbol<- genes_df[match(countNormalisedtable2$geneID, genes_df$Geneid),]$symbol
countNormalised2$symbol<-genes_df[match(countNormalised2$gene, genes_df$Geneid),]$symbol

countNormalisedtable2<-countNormalisedtable2 %>% mutate(ID=paste0(symbol,celltype,TimeCondtion,age))
countNormalisedtable<-countNormalisedtable %>% mutate(ID=paste0(symbol,celltype,TimeCondtion,age))
sortie_merge<-merge(countNormalisedtable %>% dplyr::select(ID,day,meancountnormalised),sortie_full, by="ID")
sortie_merge2<-merge(countNormalisedtable2 %>% dplyr::select(ID,meancountnormalised_t,meancountnormalised_t.1),sortie_full, by="ID")
saveRDS(sortie_merge, paste0(dynamicD,"Young/", DEG_table, "Merge_dynamic_full.rds")) 
saveRDS(sortie_merge, paste0(dynamicD,"Old/", DEG_table, "Merge_dynamic_full.rds")) 
saveRDS(sortie_merge2, paste0(dynamicD,"Young/", DEG_table, "Merge_dynamic_full2.rds")) 
saveRDS(sortie_merge2, paste0(dynamicD,"Old/", DEG_table, "Merge_dynamic_full2.rds")) 


if (volcanoneeded){
  for ( A in unique(metadata$age)) {
    sortie_full <- read.table(paste0(dynamicD,A,"/",DEG_table, A,"_dynamic_full.csv") ,sep="\t", header=T)
    names(colorsTime) = sort(unique(sortie_full$TimeCondtion))
    sortie_full <- sortie_full %>% mutate (type = factor( sortie_full$type , levels =orderTypecelldyn ))
    
    # set aesthetics data:
    sortie_full <- sortie_full %>% mutate(DEsignificant=case_when(
      padj <= 0.05 & log2FoldChange >= 1.2 ~ "Signi UP" ,
      padj <= 0.05 & log2FoldChange <= -1.2 ~ "Signi DOWN",
      TRUE  ~ "Not Signi"
    ))
    sortie_full
    textdata<-lapply(unique(sortie_full$type), function(t) lapply(unique(sortie_full %>% filter(type==t) %>% select(TimeCondtion) %>% unlist()), function(c)  rbind(sortie_full %>% filter(type==t,TimeCondtion==c) %>% filter(padj<=0.005) %>% arrange(log2FoldChange) %>% head(n=5),sortie_full %>% filter(type==t,TimeCondtion==c) %>% filter(padj<=0.005) %>% arrange(log2FoldChange) %>% tail(n=5)) ) )
    textdata<- bind_rows(textdata)
    
    g <-ggplot(sortie_full, aes(x=log2FoldChange, y = -log10(padj+1e-20),color=DEsignificant)) +
      geom_point(aes(color=DEsignificant),size=.3) +
      scale_color_manual(values=c("lightgrey",brewer.pal(10, "RdBu")[9],brewer.pal(10, "RdBu")[2] )) +
      geom_vline(xintercept = c(1.2,-1.2), data=,color= "black",
                 linetype="dashed", size=.2) +
      geom_hline(yintercept = -log10(0.05+1e-20), data=,color= "black",
                 linetype="dashed", size=.2) +
      facet_grid_blank(vars(type),vars(TimeCondtion), drop = FALSE) + 
      theme_calc() +
      theme(panel.grid.major=element_blank()) +
      geom_text_repel(
        data= textdata, aes(label=symbol, fill=DEsignificant), size=2) +
      labs(title=paste0("Dynamic (Day(t) vs Day(t+1)) across type in ",A), subtitle = "Signi = Differential gene padj < 0.05, -1.2<log2FoldChange>1.2",
           caption="vertical lines: ABS(log2FoldChange)=1.2 
            labels on top 5 log2FoldChange UP and Down for genes padj < 0.005")
    g2 <- ggplot_gtable(ggplot_build(g))
    stripRowName <- which(grepl('strip-', g2$layout$name))
    k <- 1
    fills <- c(colorsTimeCondition,colorsTypedyn)
    for (i in stripRowName) {
      j <- which(grepl('rect', g2$grobs[[i]]$grobs[[1]]$childrenOrder))
      g2$grobs[[i]]$grobs[[1]]$children[[j]]$gp$fill <- fills[k]
      k <- k+1
    }
    png(paste0(dynamicD,A,"/PlotsDEG/volcano_dynamic_",A,".png"),units = "in", width=14, height= 7, res = 300, family = "Arial")
    grid::grid.draw(g2)
    dev.off()
    tiff(paste0(dynamicD,A,"/PlotsDEG/volcano_dynamic_",A,".tiff"),units = "in", width=14, height= 7, res = 300, family = "Arial")
    grid::grid.draw(g2)
    dev.off()
    pdf(paste0(dynamicD,A,"/PlotsDEG/volcano_dynamic_",A,".pdf"), width=14, height=10)
    grid::grid.draw(g2)
    dev.off()
    dev.off()
    
  }
}

#####
# DE analysis
# dynamic analyse old vs young
# by cell-type  t-(t-1)
##### 
countNormalised3<-list()
countNormalisedtable_replicat<-data.frame(geneID = character(),age=character(),celltype=character(),TimeCondtion=character(),day=character(), replicat=character(),meancountnormalised=character())
countNormalisedtable<-data.frame(geneID = character(),age=character(),celltype=character(),TimeCondtion=character(),day=character(),meancountnormalised=character())
countNormalisedtable2<-data.frame(geneID = character(),age=character(),celltype=character(),TimeCondtion=character(),meancountnormalised_t=character(),"meancountnormalised_t+1"=character())

countzscore3<-list()
countzscoretable_replicat<-data.frame(geneID = character(),age=character(),celltype=character(),TimeCondtion=character(),day=character(),replicat=character(),zscore_count=character())
countzscoretable<-data.frame(geneID = character(),age=character(),celltype=character(),TimeCondtion=character(),day=character(),zscore_count=character())
countzscoretable2<-data.frame(geneID = character(),age=character(),celltype=character(),TimeCondtion=character(),zscore_count_t=character(),"zscore_count_t+1"=character())


timeline<-c("D0","D2","D4","D7")
typesv<-lapply( unique(metadata$day), function(d) unique(metadata[metadata$day == d,]$type) )
names(typesv)<-unique(metadata$day)
sortie_full <- data.frame("baseMean"= double(), "log2FoldChange"= double(), 
                          "lfcSE"= double(),  "stat" = double(), 
                          "pvalue" = double(), "padj" = double(),
                          "id" = character(),"TimeCondtion" = character(), "type"=character())
out_softfilt<-data.frame()
for ( t in 1:(length(timeline)-1)){
  print(t)
  D<-timeline[t]
  D2<-timeline[t+1]
  typeShared<-intersect(typesv[[D]], typesv[[D2]])
  fmat_D<- as.data.frame(fmat) %>% dplyr::select(contains(c(D,D2))) %>% dplyr::select(contains(typeShared))  %>% as.matrix()
  metadata_D <- metadata %>% dplyr::filter( day==D|day==D2, type %in% typeShared)    %>% mutate(TimeCondition=ifelse(day == D , "TimeRef", "TimePlusOne")) %>% arrange(factor(newname,levels = colnames(fmat_D)))
  
  deseqdataset.dynamic <- DESeqDataSetFromMatrix(countData = fmat_D,
                                                 colData = metadata_D,
                                                 design= ~  age +TimeCondition + age:TimeCondition)
  
  # set "Young" as reference
  colData(deseqdataset.dynamic)$TimeCondtion <- factor(colData(deseqdataset.dynamic)$TimeCondition,
                                                       levels=c("TimeRef", "TimePlusOne"))
  colData(deseqdataset.dynamic)$age <- factor(colData(deseqdataset.dynamic)$age,
                                              levels=c("Young", "Old"))
  
  
  
  for (t in typeShared){
    print(paste0(D,' vs ',D2,' in ', t))
    # do subsets as needed + deseq analysis
    deseqdataset.dynamic.daytype <- deseqdataset.dynamic[,deseqdataset.dynamic$type == t ]
    deseqdataset.dynamic.daytype <- DESeq(deseqdataset.dynamic.daytype, test="Wald", full= ~ age +TimeCondition + age:TimeCondition)
    
    countNormalised<-DESeq2::counts(deseqdataset.dynamic.daytype, normalize =T)
    countNormalised<-as.data.frame(countNormalised)
    countNormalised$id<-countNormalised %>% rownames()
    nbreplicat_DY<-length(colnames(as.data.frame(countNormalised) %>% dplyr::select(contains(paste0("Young.",t,".",D)))))
    nbreplicat_D2Y<-length(colnames(as.data.frame(countNormalised) %>% dplyr::select(contains(paste0("Young.",t,".",D2)))))
    nbreplicat_DO<-length(colnames(as.data.frame(countNormalised) %>% dplyr::select(contains(paste0("Old.",t,".",D)))))
    nbreplicat_D2O<-length(colnames(as.data.frame(countNormalised) %>% dplyr::select(contains(paste0("Old.",t,".",D2)))))
    
    total_replicat<-nbreplicat_DY+nbreplicat_D2Y+nbreplicat_DO+nbreplicat_D2O
    nbreplicatO<-nbreplicat_DO+nbreplicat_D2O
    nbreplicatY<-nbreplicat_DY+nbreplicat_D2Y
    
    countNormalisedtable_replicat_tempo<-data.frame(geneID=c(rep(countNormalised$id,total_replicat)),
                                                    age=c(rep("Young", dim(countNormalised)[1]*nbreplicatY),rep("Old", dim(countNormalised)[1]*nbreplicatO)),
                                                    celltype=rep(t, dim(countNormalised)[1]*total_replicat),
                                                    TimeCondtion=rep(paste0(D,'_',D2), dim(countNormalised)[1]*total_replicat),
                                                    day=c(rep(D, dim(countNormalised)[1]*nbreplicat_DY),rep(D2, dim(countNormalised)[1]*nbreplicat_D2Y),rep(D, dim(countNormalised)[1]*nbreplicat_DO),rep(D2, dim(countNormalised)[1]*nbreplicat_D2O)),
                                                    replicat=c(rep(c(1:nbreplicat_DY,1:nbreplicat_D2Y,1:nbreplicat_DO,1:nbreplicat_D2O),dim(countNormalised)[1])),
                                                    countnormalised=c(lapply(1:nbreplicat_DY, function(R) countNormalised[,paste0("Young.",t,".",D,"_",R)]) %>% unlist(), 
                                                                      lapply(1:nbreplicat_D2Y, function(R) countNormalised[,paste0("Young.",t,".",D2,"_",R)]) %>% unlist(), 
                                                                      lapply(1:nbreplicat_DO, function(R) countNormalised[,paste0("Old.",t,".",D,"_",R)]) %>% unlist(), 
                                                                      lapply(1:nbreplicat_D2O, function(R) countNormalised[,paste0("Old.",t,".",D2,"_",R)]) %>% unlist()))
    
    countNormalisedtable_tempo= data.frame(geneID=c(rep(countNormalised$id,4)),
                                           age=c(rep("Young", dim(countNormalised)[1]*2),rep("Old", dim(countNormalised)[1]*2)),
                                           celltype=rep(t, dim(countNormalised)[1]*4),
                                           TimeCondtion=rep(paste0(D,'_',D2), dim(countNormalised)[1]*4),
                                           day=c(rep(D, dim(countNormalised)[1]),rep(D2, dim(countNormalised)[1]),rep(D, dim(countNormalised)[1]),rep(D2, dim(countNormalised)[1])),
                                           meancountnormalised=c(rowMeans(countNormalised[,str_detect(colnames(countNormalised),paste0("Young.",t,".",D))]),rowMeans(countNormalised[,str_detect(colnames(countNormalised),paste0("Young.",t,".",D2))]),rowMeans(countNormalised[,str_detect(colnames(countNormalised),paste0("Old.",t,".",D))]),rowMeans(countNormalised[,str_detect(colnames(countNormalised),paste0("Old.",t,".",D2))])))
    countNormalisedtable_tempo2= data.frame(geneID=rep(rep(countNormalised$id,2)),
                                            age=c(rep("Young", dim(countNormalised)[1]),rep("Old", dim(countNormalised)[1])),
                                            celltype=rep(t, dim(countNormalised)[1]*2),
                                            TimeCondtion=rep(paste0(D,'_',D2), dim(countNormalised)[1]*2),
                                            meancountnormalised_t=c(rowMeans(countNormalised[,str_detect(colnames(countNormalised),paste0("Young.",t,".",D))]),rowMeans(countNormalised[,str_detect(colnames(countNormalised),paste0("Old.",t,".",D))])),
                                            "meancountnormalised_t+1"=c(rowMeans(countNormalised[,str_detect(colnames(countNormalised),paste0("Young.",t,".",D2))]),rowMeans(countNormalised[,str_detect(colnames(countNormalised),paste0("Old.",t,".",D2))])))
    
    countNormalisedtable_replicat<-rbind(countNormalisedtable_replicat,countNormalisedtable_replicat_tempo)
    countNormalisedtable<- rbind(countNormalisedtable,countNormalisedtable_tempo)
    countNormalisedtable2<- rbind(countNormalisedtable2,countNormalisedtable_tempo2)
    
    
    zscore <- t(scale(t(assay(vst(deseqdataset.dynamic.daytype))))) %>% data.frame() %>% rownames_to_column(var="id")
    zscore$symbol<-genes_df[match(zscore$id, genes_df$Geneid),]$symbol
    t_2=str_replace(t,"-",".")
    
    countzscoretable_replicat_tempo<-data.frame(geneID=c(rep(zscore$id,total_replicat)),
                                                age=c(rep("Young", dim(zscore)[1]*nbreplicatY),rep("Old", dim(zscore)[1]*nbreplicatO)),
                                                celltype=rep(t, dim(zscore)[1]*total_replicat),
                                                TimeCondtion=rep(paste0(D,'_',D2), dim(zscore)[1]*total_replicat),
                                                day=c(rep(D, dim(zscore)[1]*nbreplicat_DY),rep(D2, dim(zscore)[1]*nbreplicat_D2Y),rep(D, dim(zscore)[1]*nbreplicat_DO),rep(D2, dim(zscore)[1]*nbreplicat_D2O)),
                                                replicat=c(rep(c(1:nbreplicat_DY,1:nbreplicat_D2Y,1:nbreplicat_DO,1:nbreplicat_D2O),dim(zscore)[1])),
                                                countnormalised=c(lapply(1:nbreplicat_DY, function(R) zscore[,paste0("Young.",t_2,".",D,"_",R)]) %>% unlist(), 
                                                                  lapply(1:nbreplicat_D2Y, function(R) zscore[,paste0("Young.",t_2,".",D2,"_",R)]) %>% unlist(), 
                                                                  lapply(1:nbreplicat_DO, function(R) zscore[,paste0("Old.",t_2,".",D,"_",R)]) %>% unlist(), 
                                                                  lapply(1:nbreplicat_D2O, function(R) zscore[,paste0("Old.",t_2,".",D2,"_",R)]) %>% unlist()))
    
    zcore_tempo= data.frame(geneID=c(rep(zscore$id,4)),
                            age=c(rep("Young", dim(zscore)[1]*2),rep("Old", dim(zscore)[1]*2)),
                            celltype=rep(t, dim(zscore)[1]*4),
                            TimeCondtion=rep(paste0(D,'_',D2), dim(zscore)[1]*4),
                            day=c(rep(D, dim(zscore)[1]),rep(D2, dim(zscore)[1]),rep(D, dim(zscore)[1]),rep(D2, dim(zscore)[1])),
                            zscore_count=c(rowMeans(zscore[,str_detect(colnames(zscore),paste0("Young.",t_2,".",D))]),rowMeans(zscore[,str_detect(colnames(zscore),paste0("Young.",t_2,".",D2))]),rowMeans(zscore[,str_detect(colnames(zscore),paste0("Old.",t_2,".",D))]),rowMeans(zscore[,str_detect(colnames(zscore),paste0("Old.",t_2,".",D2))])))
    zcore_tempo2= data.frame(geneID=rep(rep(zscore$id,2)),
                             age=c(rep("Young", dim(zscore)[1]),rep("Old", dim(zscore)[1])),
                             celltype=rep(t, dim(zscore)[1]*2),
                             TimeCondtion=rep(paste0(D,'_',D2), dim(zscore)[1]*2),
                             zscore_count_t=c(rowMeans(zscore[,str_detect(colnames(zscore),paste0("Young.",t_2,".",D))]),rowMeans(zscore[,str_detect(colnames(zscore),paste0("Old.",t_2,".",D))])),
                             "zscore_count_t+1"=c(rowMeans(zscore[,str_detect(colnames(zscore),paste0("Young.",t_2,".",D2))]),rowMeans(zscore[,str_detect(colnames(zscore),paste0("Old.",t_2,".",D2))])))
    countzscoretable_replicat<-rbind(countzscoretable_replicat,countzscoretable_replicat_tempo)
    countzscoretable<- rbind(countzscoretable,zcore_tempo)
    countzscoretable2<- rbind(countzscoretable2,zcore_tempo2)
    
    # Add Shrink log2 fold changes
    deres <- lfcShrink(deseqdataset.dynamic.daytype, coef="age_Old_vs_Young", type="apeglm", parallel=T)
    deres$id <- rownames(deres)
    deres$TimeCondtion = paste0(D,"_",D2)
    deres$type = t
    
    sortie_full <- rbind(sortie_full, as_tibble(deres))
    countNormalised3[[paste0(D,'_',D2,'_', t)]]<-merge(zscore,as.data.frame(deres), by= "id" )
    
    deres_f <- as_tibble(deres) %>% dplyr::filter(padj <= 0.5 & abs(log2FoldChange) >= 0.1)
    out_softfilt <- rbind(out_softfilt,deres_f)
  }
}

# complete tables with gene_symbols & save table
sortie_full$symbol <- genes_df[match(sortie_full$id, genes_df$Geneid),]$symbol
sortie_full<- sortie_full %>% dplyr::filter(baseMean!=0) 
saveRDS(sortie_full, paste0(dynamicD, DEG_table, "dynamic_full.rds")) 
write.table(sortie_full, paste0(dynamicD,DEG_table, "dynamic_full.csv"), sep="\t", col.names=T, row.names = F) 
out_softfilt <- sortie_full %>% dplyr::filter(padj <= 0.5 & abs(log2FoldChange) >= 0.1)
write.table(out_softfilt, paste0(dynamicD,DEG_table, "dynamic_softfilter.csv"), sep="\t", col.names=T, row.names = F) 
out_hardfilt = out_softfilt %>% dplyr::filter(padj <= 0.05 & abs(log2FoldChange) >= 1.2)
write.table(out_hardfilt, paste0(dynamicD,DEG_table, "dynamic_hardfilter.csv"), sep="\t", col.names=T, row.names = F)



sortie_full<-sortie_full %>% mutate(ID=paste0(symbol,type,TimeCondtion))
countNormalisedtable_replicat$symbol<- genes_df[match(countNormalisedtable_replicat$geneID, genes_df$Geneid),]$symbol
countNormalisedtable$symbol<- genes_df[match(countNormalisedtable$geneID, genes_df$Geneid),]$symbol
countNormalisedtable2$symbol<- genes_df[match(countNormalisedtable2$geneID, genes_df$Geneid),]$symbol
countNormalisedtable2<-countNormalisedtable2 %>% mutate(ID=paste0(symbol,celltype,TimeCondtion))
countNormalisedtable<-countNormalisedtable %>% mutate(ID=paste0(symbol,celltype,TimeCondtion))
countNormalisedtable_replicat<- countNormalisedtable_replicat %>% mutate(ID=paste0(symbol,celltype,TimeCondtion))

countNormalisedtable_replicat
sortie_merge<-merge(countNormalisedtable %>% dplyr::select(ID,day,meancountnormalised,age),sortie_full, by="ID")
sortie_merge2<-merge(countNormalisedtable2 %>% dplyr::select(ID,age,meancountnormalised_t,meancountnormalised_t.1),sortie_full, by="ID")
saveRDS(sortie_merge, paste0(dynamicD, DEG_table, "Merge_dynamic_full.rds")) 
saveRDS(sortie_merge2, paste0(dynamicD, DEG_table, "Merge_dynamic_full2.rds")) 

if (volcanoneeded){
  names(colorsTime) = sort(unique(sortie_full$TimeCondtion))
  sortie_full <- readRDS(paste0(dynamicD, DEG_table, "dynamic_full.rds") )
  sortie_full <- sortie_full %>% mutate (type = factor( sortie_full$type , levels =orderTypecelldyn ))
  
  # set aesthetics data:
  sortie_full <- sortie_full %>% mutate(DEsignificant=case_when(
    padj <= 0.05 & log2FoldChange >= 1.2 ~ "Signi UP" ,
    padj <= 0.05 & log2FoldChange <= -1.2 ~ "Signi DOWN",
    TRUE  ~ "Not Signi"
  ))
  
  textdata<-lapply(unique(sortie_full$type), function(t) lapply(unique(sortie_full %>% filter(type==t) %>% select(TimeCondtion) %>% unlist()), function(c)  rbind(sortie_full %>% filter(type==t,TimeCondtion==c) %>% filter(padj<=0.005) %>% arrange(log2FoldChange) %>% head(n=5),sortie_full %>% filter(type==t,TimeCondtion==c) %>% filter(padj<=0.005) %>% arrange(log2FoldChange) %>% tail(n=5)) ) )
  textdata<- bind_rows(textdata)
  
  g <-ggplot(sortie_full, aes(x=log2FoldChange, y = -log10(padj+1e-20),color=DEsignificant)) +
    geom_point(aes(color=DEsignificant),size=.3) +
    scale_color_manual(values=c("lightgrey",brewer.pal(10, "RdBu")[9],brewer.pal(10, "RdBu")[2] )) +
    geom_vline(xintercept = c(1.2,-1.2), data=,color= "black",
               linetype="dashed", size=.2) +
    facet_grid_blank(vars(type),vars(TimeCondtion), drop = FALSE) + 
    theme_calc() +
    theme(panel.grid.major=element_blank()) +
    
    geom_text_repel(
      data= textdata,
      aes(label=symbol, fill=DEsignificant),
      size=2,
      segment.size = .1,
      force=2, force_pull = 2,
      max.overlaps=15
    ) +
    labs(title="Old vs Young dynamic (Day(t) vs Day(t+1)) across type", subtitle = "Signi = Differential gene padj < 0.05, -1.2<log2FoldChange>1.2",
         caption="vertical lines: ABS(log2FoldChange)=1.2 
            labels on top 5 log2FoldChange UP and Down for genes padj < 0.005")
  
  g2 <- ggplot_gtable(ggplot_build(g))
  stripRowName <- which(grepl('strip-', g2$layout$name))
  k <- 1
  fills <- c(colorsTimeCondition,colorsTypedyn)
  for (i in stripRowName) {
    j <- which(grepl('rect', g2$grobs[[i]]$grobs[[1]]$childrenOrder))
    g2$grobs[[i]]$grobs[[1]]$children[[j]]$gp$fill <- fills[k]
    k <- k+1
  }
  png(paste0(dynamicD,"PlotsDEG/volcano_dynamic.png"),units = "in", width=14, height= 7, res = 300, family = "Arial")
  grid::grid.draw(g2)
  dev.off()
  tiff(paste0(dynamicD,"PlotsDEG/volcano_dynamic.tiff"),units = "in", width=14, height= 7, res = 300, family = "Arial")
  grid::grid.draw(g2)
  dev.off()
  pdf(paste0(dynamicD,"PlotsDEG/volcano_dynamic.pdf"), width=14, height=10)
  grid::grid.draw(g2)
  dev.off()
  
}
