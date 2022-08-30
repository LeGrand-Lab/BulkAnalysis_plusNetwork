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
#batchinfo<-read.csv("data/batchesinfo.csv", sep = "\t")
#metadata<-merge(metadata,batchinfo)
#rownames(metadata)<-metadata$newname
fullDEsta = readRDS(paste0(odir, "TableDEG/DEG_table_static_full.rds"))
signiDEgene <- fullDEsta %>% dplyr::filter(padj<=0.05)

tableGSEA<-"GSEA/TableGSEA/"
fullGSEAconcatfile<-paste0(odir,tableGSEA,"GSEA_table_static_full.csv" )
fullGSEAconcat<-read.csv(fullGSEAconcatfile,sep = " " )
GSEAsigni<-fullGSEAconcat %>% dplyr::filter(padj<=0.05) %>% mutate(sens=ifelse(NES >0,"Up","Down"))
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
#                     ,Batch = factor(metadata$batch,levels = unique(metadata$batch)))
rownames(myannot) = rownames(metadata)

colorsTime = c(brewer.pal(9,"BuPu")[3],brewer.pal(9,"BuPu")[5],brewer.pal(9,"BuPu")[7],brewer.pal(9,"BuPu")[9])
names(colorsTime) = sort(unique(myannot$Time))
colorsAge = c("#ffc35d","#019190")
names(colorsAge) = c("Young","Old")
colorbatch<-c("white","lightgrey","darkgrey","black")
names(colorbatch)<-unique(myannot$Batch)
mycolors = list("CellType"=colorsType,"Time" = colorsTime, "Age" =colorsAge  )
                #,"Batch"=colorbatch)

##
# Correlation matrice , spearman test
##
# =========================================================================
cor.mat <- cor(fTPM, method="spearman")



heatmap_lib<-pheatmap::pheatmap(cor.mat, fontsize = 8,
    color = colorRampPalette(c("white","orange"))(100),
     cluster_rows = TRUE,
     cluster_cols = TRUE, clustering_distance_rows = "euclidean",
     clustering_distance_cols = "euclidean",  annotation_col = myannot, annotation_row =myannot,
     annotation_colors=mycolors, fontsize_row = 5, fontsize_col = 5,
     show_rownames = F, show_colnames = F,legend_labels = "1 - rho (Spearman correlation)",
     main= "Distance matrix Spearman correlation")
heatmap_lib
save_plot(paste0(ofig,"Heatmap_libraries.png"),heatmap_lib, base_width=12, base_height= 8)


lgdType<-Legend(at=levels(myannot$CellType),legend_gp = gpar(fill = colorsType), title="CellType",labels_gp = gpar (fontsize= 5),title_gp = gpar(fontsize = 6, fontface = "bold"))
lgdAge<-Legend(at=levels(myannot$Age),legend_gp = gpar(fill = colorsAge), title="Age",labels_gp = gpar (fontsize= 5),title_gp = gpar(fontsize = 6, fontface = "bold"))
lgdTime<-Legend(at=levels(myannot$Time),legend_gp = gpar(fill = colorsTime ), title="Time",labels_gp = gpar (fontsize= 5),title_gp = gpar(fontsize = 6, fontface = "bold"))
pd = packLegend(lgdType, lgdTime, lgdAge, direction = "horizontal")
pd2 = packLegend(lgdType, lgdTime,lgdAge,direction = "vertical")


ha_top<-HeatmapAnnotation(  Age = myannot$Age,Time =myannot$Time,Celltype = myannot$CellType, col = list( Celltype =colorsType, Age = colorsAge, Time = colorsTime), simple_anno_size = unit(2,"mm"), show_legend = F, annotation_name_gp = gpar(fontsize = 5) )
ha_row<-rowAnnotation(Age = myannot$Age,Time =myannot$Time,Celltype = myannot$CellType, col = list( Celltype =colorsType, Age = colorsAge, Time = colorsTime),  simple_anno_size = unit(2,"mm"),show_legend=F,annotation_name_gp = gpar(fontsize = 5) )
heatmap_cor<-ComplexHeatmap::Heatmap( cor.mat, top_annotation=ha_top, left_annotation = ha_row ,   colorRampPalette(c("white","orange"))(100) , show_row_names = F, show_column_names = F,
                         heatmap_legend_param = list( title = "1 - rho (Spearman correlation)", direction="horizontal",title_position="topcenter",     title_gp = gpar(fontsize = 6, fontface = "bold"), labels_gp = gpar(fontsize = 5))  , column_title ="Distance matrix Spearman correlation",column_title_gp = gpar(fontsize = 10, fontface = "bold"))

png(paste0(ofig,"Heatmap_libraries1.png"),units = "in", width=10, height= 4, res = 300, family = "Arial")
ComplexHeatmap::draw(heatmap_cor,heatmap_legend_side="bottom", annotation_legend_list=pd) 
dev.off()
png(paste0(ofig,"Heatmap_libraries2.png"),units = "in", width=4, height= 5.5, res = 300, family = "Arial")
ComplexHeatmap::draw(heatmap_cor,heatmap_legend_side="bottom", annotation_legend_list=pd2) 
dev.off()


##
# II. **** PCA  ****
##
# For PCA: this uses vst ,  requires DESeq2 object
# =========================================================================
#PCA on gobal table
fmat <- readRDS("data/prefiltered_counts.rds")
metadata <- readRDS(metadata.rds)

ds.o <- DESeqDataSetFromMatrix(countData=fmat,
                                colData=metadata,
                                design= ~type + time + age:time)
dim(counts(ds.o))
vst <- varianceStabilizingTransformation(ds.o)
pca_all <- plotPCA(vst, intgroup=c("type", "time","age"), returnData=T)
percentVar.subset <- round(100* attr(pca_all,"percentVar"))
pca_all.perc.P1 <- percentVar.subset[1]
pca_all.perc.P2 <- percentVar.subset[2]

# PCA on each cell cell type
pca_df_bytype <- data.frame("PC1"=numeric(),  "PC2"=numeric(), "group"=character(), "time"=character(),"age" =character(),"name"=character() ,"type2"=character())
perc.PC1.type <- c()
perc.PC2.type <- c() 
for (k in unique(myannot$CellType) ){
  subset <- ds.o[,ds.o$type==k]
  subset.vst <- varianceStabilizingTransformation(subset)
  tmp <- plotPCA(subset.vst, intgroup=c("time","age"), returnData=T)
  tmp$type2 <- k
  pca_df_bytype <- rbind(pca_df_bytype,tmp)
  percentVar.subset <- round(100* attr(tmp,"percentVar"))
  perc.PC1.type <- c(perc.PC1.type, percentVar.subset[1])
  perc.PC2.type <- c(perc.PC2.type, percentVar.subset[2])
}

pca_df_bytype$time<-factor(pca_df_bytype$time,levels =  levels(myannot$Time))
pca_df_bytype$type2<-factor(pca_df_bytype$type2 ,levels =  levels(myannot$CellType))
#PCA on each time

pca_df_bytime <- data.frame("PC1"=numeric(), "PC2"=numeric(),"group"=character(), "time"=character(), "age" =character(),"name"=character())
perc_PC1_time <- c()
perc_PC2_time <- c()
days <- c("D0","D2","D4","D7")
for (day in days) {
  subset <- ds.o[,ds.o$time==day]
  subset.vst <- varianceStabilizingTransformation(subset)
  tmp <- plotPCA(subset.vst, intgroup=c("age","type"), returnData=T)
  tmp <- cbind(tmp,rep(day,dim(tmp)[1]))
  pca_df_bytime <- rbind(pca_df_bytime,tmp)
  percentVar.subset <- round(100* attr(tmp,"percentVar"))
  perc_PC1_time <- c(perc_PC1_time, percentVar.subset[1])
  perc_PC2_time <- c(perc_PC2_time, percentVar.subset[2])
}
pca_df_bytime<-pca_df_bytime %>% mutate(time=`rep(day, dim(tmp)[1])`)
pca_df_bytime$time<-factor(pca_df_bytime$time,levels =  levels(myannot$Time))
pca_df_bytime$type<-factor(pca_df_bytime$type ,levels =  levels(myannot$CellType))

# PCA on each time and cell type 
metadata_daytype  <- metadata %>% mutate(timetype=paste0(time,".",type))
pca_df_daytype <- data.frame("PC1"=numeric(), "PC2"=numeric(),"group"=character(), "time"=character(), "age" =character(),"name"=character() ,"timetype"=character())
perc.PC1_daytype <- c()
perc.PC2.daytype <- c()
tp <- unique(metadata_daytype$timetype)
for (k in tp){
  tup <- str_split(k,"\\.")[[1]]
  subset <- ds.o[,ds.o$type==tup[2] & ds.o$time==tup[1]]
  subset.vst <- varianceStabilizingTransformation(subset)
  tmp <- plotPCA(subset.vst, intgroup=c("age"), returnData=T)
  tmp$timetype <- k
  pca_df_daytype <- rbind(pca_df_daytype,tmp)
  percentVar.subset <- round(100* attr(tmp,"percentVar"))
  perc.PC1_daytype <- c(perc.PC1_daytype, percentVar.subset[1])
  perc.PC2.daytype <- c(perc.PC2.daytype, percentVar.subset[2])
}

pca_df_daytype$time <- unname(sapply(pca_df_daytype$timetype, 
                              function(x){ str_split(x,"\\.")[[1]][1]}))
pca_df_daytype$type <- unname(sapply(pca_df_daytype$timetype, 
                              function(x){ str_split(x,"\\.")[[1]][2]}))
pca_df_daytype$type <- factor(pca_df_daytype$type, levels = levels(myannot$CellType))
pca_df_daytype$age<-factor(pca_df_daytype$age, levels = levels(myannot$Age))

## Productio plot

margin(t = 0, r = 0, b = 0, l = 0, unit = "cm")
pca_all<-pca_all %>% mutate(type=factor(pca_all$type,level=orderTypecell))

facet1 <- ggplot(pca_all, aes(x=PC1, y=PC2,  color=type)) +
  geom_point(size=0.8, show.legend = F) + 
  stat_ellipse(aes(group=type), size=0.2, show.legend = F) + 
  scale_color_manual(values=colorsType) +
  theme_light() +
  theme(axis.title = element_text(size=5, face = "bold"),axis.text=element_text(size=3),legend.title=element_text(size=6, face = "bold"),legend.text=element_text(size=4),title=element_text(size=7, face = "bold"),
        legend.spacing= unit(0.1,"points"),legend.spacing.y=unit(0.01,"mm"),legend.key.size = unit(3, "mm"),panel.grid =element_line(color="white"))+
  xlab(paste0("PC1: ", pca_all.perc.P1,"% variance" )) +
  ylab(paste0("PC2: ", pca_all.perc.P2,"% variance" )) +
  ggtitle("The first two components of the PCA on whole dataset")
facet1

lgdType_grob=grid.grabExpr(draw(lgdType)) 
facet1_ld<-grid.arrange(facet1,lgdType_grob, ncol=2,widths=c(4,1))
facet1_ld

save_plot(paste0(ofig,"PCA_grouping_by_age_time_typecell.png"),facet1_ld)

facets2 <- ggplot(pca_df_bytype, aes(x=PC1, y=PC2, color=time,fill=time, shape=age)) +
  geom_point(size=0.8,aes(shape=age), color="grey") + 
  scale_fill_manual(values=colorsTime) +
  scale_color_manual(values=colorsTime) +
  scale_shape_manual(values=c(21,25))+
  stat_ellipse(aes(group=time),size=0.2) + 
  theme_light() +
  theme(axis.title = element_text(size=5, face = "bold"),axis.text=element_text(size=3),legend.title=element_text(size=6, face = "bold"),legend.text=element_text(size=4),title=element_text(size=7, face = "bold"),
        legend.spacing= unit(0.1,"points"),legend.spacing.y=unit(0.01,"points"),strip.text.x = element_text(size = 4, face = "bold", color="white"),strip.text.y = element_text(size=4, face = "bold", color="black"),legend.key.size = unit(3, "mm"),panel.spacing = unit(0.5,"mm"),panel.grid =element_line(color="white"))+
  facet_grid(vars(type2))

facets2
g2 <- ggplot_gtable(ggplot_build(facets2))
stripRowName <- which(grepl('strip-r', g2$layout$name))
k <- 1
fills <- c(colorsType)
for (i in stripRowName) {
  j <- which(grepl('rect', g2$grobs[[i]]$grobs[[1]]$childrenOrder))
  g2$grobs[[i]]$grobs[[1]]$children[[j]]$gp$fill <- fills[k]
  k <- k+1
}
grid::grid.draw(g2)

save_plot(paste0(ofig,"PCA_typecell_grouping_age_time.png"),g2, base_width=8, base_height= 3)

facets2.5 <- ggplot(pca_df_bytype, aes(x=PC1, y=PC2, color=time)) +
  geom_point(size=0.8, show.legend = F) + 
  scale_color_manual(values=colorsTime) +
  scale_x_continuous(n.breaks = 3)+
  stat_ellipse(aes(group=time),size=0.2, show.legend = F) + 
  theme_bw() +
  theme(axis.title = element_text(size=5, face = "bold"),axis.text=element_text(size=3),legend.title=element_text(size=6, face = "bold"),legend.text=element_text(size=4),title=element_text(size=7, face = "bold"),
        legend.spacing= unit(0.1,"points"),legend.spacing.y=unit(0.01,"points"),strip.text.x = element_text(size = 4, face = "bold", color="black"),strip.text.y = element_text(size=4, face = "bold", color="black"),legend.key.size = unit(3, "mm"),panel.spacing = unit(0.5,"mm"),panel.grid =element_line(color="white"))+
  xlab("PC1") +
  ylab("PC2") +
  facet_grid(. ~ type2)+
  ggtitle("PCA on sub-matrices at each type")


g2.5 <- ggplot_gtable(ggplot_build(facets2.5))
stripRowName <- which(grepl('strip-t', g2.5$layout$name))
k <- 1
fills <- c(colorsType)
for (i in stripRowName) {
  j <- which(grepl('rect', g2.5$grobs[[i]]$grobs[[1]]$childrenOrder))
  g2.5$grobs[[i]]$grobs[[1]]$children[[j]]$gp$fill <- fills[k]
  k <- k+1
}
grid::grid.draw(g2.5)

lgdTime_grob=grid.grabExpr(draw(lgdTime)) 
facet2.5_ld<-grid.arrange(g2.5,lgdTime_grob, ncol=2,widths=c(5,0.6))
facet2.5_ld

pdf(paste0(ofig,"pca_bytype.pdf"), height=4, width=12)
facets2.5
dev.off()

# ==== time splitted:
facet3<- ggplot(pca_df_bytime, aes(x=PC1,y=PC2, color=type)) +
  geom_point(size=0.8, show.legend = F) + 
  scale_color_manual(values=colorsType) +
  scale_x_continuous(n.breaks = 3)+
  stat_ellipse(aes(group=type),size=0.2, show.legend = F) + 
  theme_bw() +
  theme(axis.title = element_text(size=5, face = "bold"),axis.text=element_text(size=3),legend.title=element_text(size=6, face = "bold"),legend.text=element_text(size=4),title=element_text(size=7, face = "bold"),
        legend.spacing= unit(0.1,"points"),legend.spacing.y=unit(0.01,"points"),strip.text.x = element_text(size = 4, face = "bold", color="white"),strip.text.y = element_text(size=4, face = "bold", color="black"),legend.key.size = unit(3, "mm"),panel.spacing = unit(0.5,"mm"),panel.grid =element_line(color="white"))+
  xlab("PC1") +
  ylab("PC2") +
  facet_grid(. ~ time) +
  ggtitle("PCA on vst sub-matrices at each time")
facet3

g3 <- ggplot_gtable(ggplot_build(facet3))
stripRowName <- which(grepl('strip-t', g3$layout$name))
k <- 1
fills <- c(colorsTime )
for (i in stripRowName) {
  j <- which(grepl('rect', g3$grobs[[i]]$grobs[[1]]$childrenOrder))
  g3$grobs[[i]]$grobs[[1]]$children[[j]]$gp$fill <- fills[k]
  k <- k+1
}
grid::grid.draw(g3)

facet3_ld<-grid.arrange(g3,lgdType_grob, ncol=2,widths=c(5,0.8))

pdf(paste0(ofig,"pca_bytime.pdf"))
facet3_ld
dev.off()

# ==== all (time and type) splitted

facet4<-ggplot(pca_df_daytype, aes(x=PC1,y=PC2, color=age)) +
  geom_point(size=0.2, show.legend = F) +
  scale_color_manual(values = colorsAge) +
  facet_grid_blank(vars(type),vars(time), drop = FALSE) +
  theme_bw() +
  theme(axis.title = element_text(size=5, face = "bold"),axis.text=element_text(size=3),legend.title=element_text(size=4, face = "bold"),legend.text=element_text(size=2),title=element_text(size=7, face = "bold"),
        legend.spacing= unit(0.1,"points"),legend.spacing.y=unit(0.01,"points"),strip.text.x = element_text(size = 4, face = "bold", color="white"),strip.text.y = element_text(size=3, face = "bold", color="black"),legend.key.size = unit(3, "mm"),panel.spacing = unit(0.5,"mm"),panel.grid =element_line(color="white"))+ 
  geom_mark_ellipse(expand = unit(0.6, "mm"),size=0.2, show.legend = F)+
  facet_grid_blank(vars(type),vars(time), drop = FALSE) +
  ggtitle("PCA on sub-matrices at each time and type cell")
facet4

g4 <- ggplot_gtable(ggplot_build(facet4))
stripRowName <- which(grepl('strip-', g4$layout$name))
k <- 1
fills <- c(colorsTime,colorsType)
for (i in stripRowName) {
  j <- which(grepl('rect', g4$grobs[[i]]$grobs[[1]]$childrenOrder))
  g4$grobs[[i]]$grobs[[1]]$children[[j]]$gp$fill <- fills[k]
  k <- k+1
}
lgdAge_grob=grid.grabExpr(draw(lgdAge)) 
facet4_ld<-grid.arrange(g4,lgdAge_grob, ncol=2,widths=c(5,0.6))
pdf(paste0(ofig,"pca_typebytime.pdf"))
facet4_ld
dev.off()
save_plot(paste0(ofig,"pca_typebytime.png"),facet4_ld, base_width=10, base_height= 6)

png(paste0(ofig,"PCA_AllData_byCelltype_byCellTypeTime.png"),units = "in", width=10, height= 6, res = 300, family = "Arial")
facet1 / ( facets2.5 |facet4 )
dev.off()
png(paste0(ofig,"PCA_AllData_byCelltype_byCellTypeTime2.png"),units = "in", width=4, height= 5.5, res = 300, family = "Arial")
pcas_plot<-grid.arrange(facet1_ld , facet2.5_ld , facet4_ld, ncol=1,heights = c(0.25,0.25,0.5)) 
pcas_plot<-grid.arrange(facet1 , lgdType_grob , g2.5, lgdTime_grob , g4, lgdAge_grob, ncol=2, heights = c(0.25,0.25,0.5), widths=c(5,0.8)) 
dev.off()


png(paste0(ofig,"Figure1_head.png"),units = "in", width=8, height= 5.5, res = 300, family = "Arial")
CovHeatmap_grob=grid.grabExpr(draw(heatmap_cor,heatmap_legend_side="bottom", annotation_legend_list=pd2, annotation_legend_side="right",
                                   legend_grouping = "original")) 
ggpubr::ggarrange(CovHeatmap_grob,ggpubr::ggarrange(facet1_ld , facet2.5_ld , facet4_ld, ncol=1,heights = c(0.25,0.25,0.5), labels=c("B","C","D")), ncol=2, labels = "A")
dev.off()


#Violine plot of pathways enrichments

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
        panel.background = element_blank())+labs(fill = "Count pathways padj<0.05")+
  ggtitle("Distribution of pathways significant enriched Old vs Young \n according to log2Foldchange")
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
vertices <- data_edged%>% mutate(showlabel=ifelse(level == 2 , NamePath,NA)) %>% dplyr::select (Node2,size,NamePath,level,showlabel) %>% rbind(c("root1",7,"root",1,NA))
vertices[vertices$NamePath == "REPRODUCTION",]$showlabel<-""
vertices <- vertices %>% mutate(NbCondSigni=infoGSEA[match(paste0(" ", vertices$NamePath), infoGSEA$NamePath),]$NbCondSigni) %>% mutate(NbCondSigni=ifelse(is.na(NbCondSigni), 0, NbCondSigni))
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
getPalette = viridis(8, option = "mako", direction = -1 )

circles<-ggraph(mygraph, layout = 'circlepack', weight=as.numeric(NbTotPathEnrich)) + 
  geom_node_circle(aes(fill= factor(NbCondSigni,levels = 0:8) ,group = level2)) +
  geom_node_label(aes(label=showlabel,size=NbTotPathEnrich),repel = T) +
  theme_void() + ggtitle("Significant enriched pathways overview Old vs Young") +
  theme(legend.position=c(fill="bottom"),legend.title = element_text(size=6,face="bold"),legend.text=element_text(size=4), legend.key.size = unit(0.4,'cm'), plot.title=element_text(face="bold",size=10,hjust=0.5)) + 
  scale_fill_manual(name= "Number of celltype_day \n with pathway enriched",values = c("white",getPalette)) + scale_size_continuous(name="Number pathway enriched \n in this event hierarchy",range = c(1, 4))

print(circles)

png(file=paste0(ofig,"HirarchieALLpathwaysNEs2.png"),units = "in", width=10, height=8, res = 800, family = "Arial")
print(circles)
dev.off()

png(paste0(ofig,"FigurePCA_GSEA_corpus.png"),units = "in", width=10, height= 12, res = 800, family = "Arial")
ggpubr::ggarrange( ggpubr::ggarrange( CovHeatmap_grob,ggpubr::ggarrange(facet1 , g2.5 , g4, ncol=1,heights = c(0.25,0.25,0.5), labels=c("B","C","D")), ncol=2,labels=c("A",""),label.x=c(0.2,0)), ggpubr::ggarrange(circles,g5,ncol=2,labels = c("E","F"),label.x=c(0.1,0)),nrow = 2)
dev.off()

tiff(paste0(ofig,"FigurePCA_GSEA_corpus.tiff"),units = "in", width=10, height= 12, res = 800, family = "Arial")
ggpubr::ggarrange( ggpubr::ggarrange( CovHeatmap_grob,ggpubr::ggarrange(facet1 , g2.5 , g4, ncol=1,heights = c(0.25,0.25,0.5), labels=c("B","C","D")), ncol=2,labels=c("A",""),label.x=c(0.2,0)), ggpubr::ggarrange(circles,g5,ncol=2,labels = c("E","F"),label.x=c(0.1,0)),nrow = 2)
dev.off()
