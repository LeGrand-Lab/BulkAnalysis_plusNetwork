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
cor.mat <- cor(fTPM, method="spearman")

myannot = data.frame(CellType=factor(metadata$type, levels = orderTypecell),
                           Time = factor(metadata$time,levels= str_sort(unique(metadata$time))),
                     Age = factor(metadata$age, levels = str_sort(unique(metadata$age), decreasing =T)))
rownames(myannot) = rownames(metadata)

orderTypecell=c("ECs","FAPs","MuSCs","Neutrophils","Inflammatory-Mac","Resolving-Mac")
colorsType=c("#10b387ff","#3d85c6ff","#b171f1ff","#f0e442ff","#ff9900ff","#cc0000ff")
names(colorsType) = orderTypecell
colorsTime = c(brewer.pal(9,"BuPu")[3],brewer.pal(9,"BuPu")[5],brewer.pal(9,"BuPu")[7],brewer.pal(9,"BuPu")[9])
names(colorsTime) = sort(unique(myannot$Time))
colorsAge = c("#ffc35d","#019190")
names(colorsAge) = c("Young","Old")
mycolors = list("CellType"=colorsType,"Time" = colorsTime, "Age" =colorsAge )

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

png(paste0(ofig,"Figure1_corpus.png"),units = "in", width=8, height= 11, res = 300, family = "Arial")
ggpubr::ggarrange(CovHeatmap_grob,ggpubr::ggarrange(facet1_ld , facet2.5_ld , facet4_ld, ncol=1,heights = c(0.25,0.25,0.5), labels=c("B","C","D")),g5, ncol=2,nrow = 2, labels = c("A","","E"))
dev.off()
png(paste0(ofig,"Figure1_tail.png"),units = "in", width=4, height= 5.5, res = 300, family = "Arial")
grid::grid.draw(g5)
dev.off()


save_plot(paste0(ofig,"PlotDEGsigniDayTypeVioline.png"),g5,base_width = 10, base_height = 7)

### Place sankey plot
colorCombiDayDEG<-c("#fef6be","#e2f4df","#dbeaf5","#d8d9ea",
                    "#d6efa6","#c2d9ed","#ffcce4","#a0dbb5","#d5b7d7","#8895c2",
                    "#008349","#065fa3","#d42857","#7e2175","#67081d")
names(colorCombiDayDEG)<-c("D0","D2","D4","D7"
                           ,"D0.D2","D0.D4","D0.D7","D2.D4","D2.D7","D4.D7",
                           "D0.D2.D4","D0.D2.D7","D0.D4.D7","D2.D4.D7","D0.D2.D4.D7")

linksToRemove<-c("D0.D4","D0.D7","D2.D7","D0.D2.D4","D0.D2.D7","D0.D4.D7")
colorsConditions<-colorsConditions[! names(colorsConditions) %in% linksToRemove]
lgddaysDEG<-Legend(at=names(colorsConditions),legend_gp = gpar(fill = colorsConditions), title="days_DEG",labels_gp = gpar (fontsize= 5),title_gp = gpar(fontsize = 6, fontface = "bold"))
lgddaysDEG_grob=grid.grabExpr(draw(lgddaysDEG)) 
grid.arrange(lgddaysDEG_grob)

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
facet6_ld<-grid.arrange(g6,lgddaysDEG_grob,ncol=2,widths=c(7,1.5))

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
nodes<-nodes  %>% arrange(loops,type,day)
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

PantaRhei::sankey(nodes, links2, palettesCombiDayDEG,
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



