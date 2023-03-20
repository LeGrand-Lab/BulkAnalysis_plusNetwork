##
# I. TPM across libraries Pheatmap. II. RAW COUNTS --> vst --> PCA.
# johaGL 2021
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
library(scales)
library(ComplexHeatmap)
library(patchwork)
library(facetious)
ttf_import(paths = "/home/bioinfo/R/fonts")
fonts()
loadfonts(device="postscript")
setwd("~/BulkAnalysis_plusNetwork/")


ofig <- "plotsPrelim/"
#outputs:
prefil_cou <- "data/prefiltered_counts.rds"
prefil_tpm <- "data/prefiltered_TPM.rds"
metadata.rds <- "data/metadata.rds"


# PRELIMINAR :  counts and TPM matrices formatting
# =========================================================================
# load protein coding matrices :
pc.mat <- readRDS("data/protcod_counts.rds") 
pc.TPM <- readRDS("data/protcod_TPM.rds")
design <- read.table("data/design.csv",header=T,sep=",", row.names=1)

# Remove transcrits associated with gene witch associated with several genes
genes_df <- read.table("data/genesinfo.csv", sep="\t", header=T)
pc.TPM.tempo <- data.frame(pc.TPM ,symbol=genes_df[match(rownames(pc.TPM), genes_df$Geneid),]$symbol)
extractGeneDuplicate <- rownames(pc.TPM.tempo[duplicated(pc.TPM.tempo$symbol),])

pc.mat<-pc.mat[!rownames(pc.mat) %in% extractGeneDuplicate,]
pc.TPM<-pc.TPM[!rownames(pc.TPM) %in% extractGeneDuplicate,]


# make sure numeric format (as.matrix), and filter out all-zeroes rows:

fTPM <- mutate_all(pc.TPM, function(x) as.numeric(as.character(x)))
fTPM <- as.matrix(pc.TPM, dim=dim(pc.TPM))
fmat <- mutate_all(pc.mat, function(x) as.numeric(as.character(x)))
fmat <- as.matrix(pc.mat, dim=dim(pc.mat))
fTPM <- fTPM[!rowSums(fTPM)==0,]  # yes, c'est bien COUNTS criterium, OK 
fmat <- fmat[!rowSums(fmat)==0,]  
typeof(fmat[3000,50])
#[1] "integer"
badquality <- setdiff(colnames(fmat),design$sample)
# bad quality libraries (cf huong) , they were already deleted from design.csv:
#[1] "X2.422363" "X2.422364" "X2.424030" "X2.424031" "X2.428151" "X2.428100" "X2.429794" "X2.429793"
exclu <- fmat[,badquality]
fmat <- fmat[,!(colnames(fmat) %in% badquality)]
fTPM <- fTPM[,!(colnames(fTPM) %in% badquality)]
# take away REDO and All_Neg bad quality libraries as well : 
metadata <- design
metadata <- metadata[!str_detect(metadata$type,"All_Neg"),]
metadata <- metadata[!str_detect(metadata$type,"REDO"),]
#retain only good quality libraries (those included in filtered metadata):
fmat <- fmat[,colnames(fmat) %in% metadata$sample]
fTPM <- fTPM[,colnames(fTPM) %in% metadata$sample]

# rename sample names:
#  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
metadata$type <- str_replace(str_replace(metadata$type,"-RNA",""),"-REDO","")
metadata$group <- paste0(metadata$age,".",metadata$type)
newnames <- paste0(metadata$age,".", metadata$type, ".", metadata$time)
names(newnames) <- metadata$sample
newnames <- sort(newnames)
replicnbs <- tibble("newn"=newnames) %>% group_by(newn) %>% tally() #counts() deprecated
nbs <- c()  # a vector containing only numerical suffix to distinguish replicates
for (i in 1:dim(replicnbs)[1]){
  nbs <- c(nbs,seq(1:replicnbs$n[i]))
}
dfnwn <- data.frame("sample"=names(newnames),"nwn"=newnames,"replicate"=nbs)
dfnwn <- dfnwn %>% mutate(newnames=paste0(nwn,"_",replicate))
## add these new names to the design
metadata$newname = dfnwn$newnames[match(metadata$sample,dfnwn$sample)]
rownames(metadata) = metadata$newname
# and also to the matrices
precols <- colnames(fmat)
tt <- c()
for (i in 1:length(precols)){
    key <- precols[i]
    nn <- metadata[metadata$sample==key,]$newname
    if (!is.na(nn)){
      tt <- c(tt, nn) 
    }else{ tt <- c(tt, key)}
}
colnames(fmat) <- unname(tt)
colnames(fTPM) <- unname(tt)
#set same order 
metadata <- metadata[match(colnames(fmat),rownames(metadata)),] 

# end rename sample names
# save these new matrices, and metadata
saveRDS(fTPM, file=prefil_tpm)
saveRDS(fmat, file=prefil_cou)
saveRDS(metadata, file=metadata.rds)
# ==========================================================================
# END preliminar

fmat <- readRDS(prefil_cou)
fTPM <- readRDS(prefil_tpm)
metadata <- readRDS(metadata.rds)

##
# I. **** TPM libraries pheatmap ****
##
# =========================================================================
cor.mat <- cor(fTPM, method="spearman")

#batch = read.table("data/batchesinfo.csv") # to add in visuals
#metadata$batch = batch[match(metadata$sample, batch$sample),]$batch

                     
rownames(myannot) = rownames(metadata)


show_col(colorblind_pal()(8))
colorsType = c(colorblind_pal()(8)[7],colorblind_pal()(8)[8],"#f6e303",colorblind_pal()(8)[3],colorblind_pal()(8)[6],colorblind_pal()(8)[4])
names(colorsType) = unique(myannot$CellType)
colorsTime = c(brewer.pal(9,"BuPu")[3],brewer.pal(9,"BuPu")[5],brewer.pal(9,"BuPu")[7],brewer.pal(9,"BuPu")[9])
names(colorsTime) = sort(unique(myannot$Time))
colorsAge = c("#ffc35d","#019190")
names(colorsAge) = c("Young","Old")
mycolors = list("CellType"=colorsType,
              "Time" = colorsTime,
              "Age" =colorsAge )

heatmap_lib<-pheatmap::pheatmap(cor.mat, fontsize = 8,
    color = colorRampPalette(c("white","orange"))(100),
     cluster_rows = TRUE,
     cluster_cols = TRUE, clustering_distance_rows = "euclidean",
     clustering_distance_cols = "euclidean",  annotation_col = myannot, annotation_row =myannot,
     annotation_colors=mycolors, fontsize_row = 5, fontsize_col = 5,
     show_rownames = F, show_colnames = F,legend_labels = "1 - rho (Spearman correlation)",
     main= "Distance matrix Spearman correlation")
heatmap_lib
png(paste0(ofig,"Heatmap_libraries1.png"),units = "in", width=10, height= 4, res = 300, family = "Arial")
png(paste0(ofig,"Heatmap_libraries2.png"),units = "in", width=4, height= 5.5, res = 300, family = "Arial")
lgdType<-Legend(at=levels(myannot$CellType),legend_gp = gpar(fill = colorsType), title="CellType",labels_gp = gpar (fontsize= 5),title_gp = gpar(fontsize = 6, fontface = "bold"))
lgdAge<-Legend(at=levels(myannot$Age),legend_gp = gpar(fill = colorsAge), title="Age",labels_gp = gpar (fontsize= 5),title_gp = gpar(fontsize = 6, fontface = "bold"))
lgdTime<-Legend(at=levels(myannot$Time),legend_gp = gpar(fill = colorsTime ), title="Time",labels_gp = gpar (fontsize= 5),title_gp = gpar(fontsize = 6, fontface = "bold"))
pd = packLegend(lgdType, lgdTime, lgdAge, direction = "horizontal")
pd2 = packLegend(lgdType, lgdTime,lgdAge,direction = "vertical")


ha_top<-HeatmapAnnotation(  Age = myannot$Age,Time =myannot$Time,Celltype = myannot$CellType, col = list( Celltype =colorsType, Age = colorsAge, Time = colorsTime), simple_anno_size = unit(2,"mm"), show_legend = F, annotation_name_gp = gpar(fontsize = 5) )
ha_row<-rowAnnotation(Age = myannot$Age,Time =myannot$Time,Celltype = myannot$CellType, col = list( Celltype =colorsType, Age = colorsAge, Time = colorsTime),  simple_anno_size = unit(2,"mm"),show_legend=F,annotation_name_gp = gpar(fontsize = 5) )
heatmap_cor<-ComplexHeatmap::Heatmap( cor.mat, top_annotation=ha_top, left_annotation = ha_row ,   colorRampPalette(c("white","orange"))(100) , show_row_names = F, show_column_names = F,
                         heatmap_legend_param = list( title = "1 - rho (Spearman correlation)", direction="horizontal",title_position="topcenter",     title_gp = gpar(fontsize = 6, fontface = "bold"), labels_gp = gpar(fontsize = 5))  , column_title ="Distance matrix Spearman correlation",column_title_gp = gpar(fontsize = 10, fontface = "bold"))
ComplexHeatmap::draw(heatmap_cor,heatmap_legend_side="bottom", annotation_legend_list=pd) 
ComplexHeatmap::draw(heatmap_cor,heatmap_legend_side="bottom", annotation_legend_list=pd2) 
save_plot(paste0(ofig,"Heatmap_libraries.png"),heatmap_lib, base_width=12, base_height= 8)
dev.off()

# rev(brewer.pal(n = 7, name = "PRGn")))(100),

# #  testing same cor function on log10(TPM+1) data : same result
# cormatB <- cor(log10(fTPM+1))
# pheatmap(cormatB, method="spearman",fontsize = 7)
# =========================================================================
# end TPM libraries composition
 
##
# II. **** PCA  ****
##
# For PCA: this uses vst ,  requires DESeq2 object
# =========================================================================
#PCA on gobal table
ds.o <- DESeqDataSetFromMatrix(countData=fmat,
                                colData=metadata,
                                design= ~type + time + age:time)
dim(counts(ds.o))
vst <- varianceStabilizingTransformation(ds.o)
pcavst3 <- plotPCA(vst, intgroup=c("type", "time","age"), returnData=T)
percentVar.subset <- round(100* attr(pcavst3,"percentVar"))
pcavst3.perc.P1 <- percentVar.subset[1]
pcavst3.perc.P2 <- percentVar.subset[2]
# PCA on each cell cell type
pca_dfA <- data.frame("PC1"=numeric(),   
                      "PC2"=numeric(),
                      "group"=character(), "time"=character(),
                      "age" =character(),"name"=character() ,
                      "type2"=character())
perc.P1 <- c()
perc.P2 <- c() 
for (k in unique(myannot$CellType) ){
  subset <- ds.o[,ds.o$type==k]
  subset.vst <- varianceStabilizingTransformation(subset)
  tmp <- plotPCA(subset.vst, intgroup=c("time","age"), returnData=T)
  tmp$type2 <- k
  pca_dfA <- rbind(pca_dfA,tmp)
  percentVar.subset <- round(100* attr(tmp,"percentVar"))
  perc.P1 <- c(perc.P1, percentVar.subset[1])
  perc.P2 <- c(perc.P2, percentVar.subset[2])
}

pca_dfA$time<-factor(pca_dfA$time,levels =  levels(myannot$Time))
pca_dfA$type2<-factor(pca_dfA$type2 ,levels =  levels(myannot$CellType))
#PCA on each time
ds.o <- DESeqDataSetFromMatrix(countData=fmat,
                               colData=metadata,
                               design= ~age)

metadataB <- metadata
pca_dfB <- data.frame("PC1"=numeric(),   
                      "PC2"=numeric(),
                      "group"=character(), "time"=character(),
                      "age" =character(),"name"=character())
pP1 <- c()
pP2 <- c()
days <- c("D0","D2","D4","D7")
for (day in days) {
  subset <- ds.o[,ds.o$time==day]
  subset.vst <- varianceStabilizingTransformation(subset)
  tmp <- plotPCA(subset.vst, intgroup=c("age","type"), returnData=T)
  tmp <- cbind(tmp,rep(day,dim(tmp)[1]))
  pca_dfB <- rbind(pca_dfB,tmp)
  percentVar.subset <- round(100* attr(tmp,"percentVar"))
  pP1 <- c(pP1, percentVar.subset[1])
  pP2 <- c(pP2, percentVar.subset[2])
}
pca_dfB<-pca_dfB %>% mutate(time=`rep(day, dim(tmp)[1])`)
# PCA on each time and cell type 
metadataC  <- metadata %>% mutate(timetype=paste0(time,".",type))
pca_dfC <- data.frame("PC1"=numeric(),   
                      "PC2"=numeric(),
                      "group"=character(), "time"=character(),
                      "age" =character(),"name"=character() ,
                      "timetype"=character())
perc.P1 <- c()
perc.P2 <- c()
tp <- unique(metadataC$timetype)
for (k in tp){
  tup <- str_split(k,"\\.")[[1]]
  subset <- ds.o[,ds.o$type==tup[2] & ds.o$time==tup[1]]
  subset.vst <- varianceStabilizingTransformation(subset)
  tmp <- plotPCA(subset.vst, intgroup=c("age"), returnData=T)
  tmp$timetype <- k
  pca_dfC <- rbind(pca_dfC,tmp)
  percentVar.subset <- round(100* attr(tmp,"percentVar"))
  perc.P1 <- c(perc.P1, percentVar.subset[1])
  perc.P2 <- c(perc.P2, percentVar.subset[2])
}
variancesdf <- data.frame("timetype"=tp, "PC1"=perc.P1,"PC2"=perc.P2 )

pca_dfC$time <- unname(sapply(pca_dfC$timetype, 
                              function(x){ str_split(x,"\\.")[[1]][1]}))
pca_dfC$type <- unname(sapply(pca_dfC$timetype, 
                              function(x){ str_split(x,"\\.")[[1]][2]}))
pca_dfC$type <- factor(pca_dfC$type, levels = levels(myannot$CellType))
pca_dfC$age<-factor(pca_dfC$age, levels = levels(myannot$Age))
margin(t = 0, r = 0, b = 0, l = 0, unit = "cm")
facet1 <- ggplot(pcavst3, aes(x=PC1, y=PC2,  color=type)) +
  geom_point(size=0.8, show.legend = F) + 
  stat_ellipse(aes(group=type), size=0.2, show.legend = F) + 
  scale_color_manual(values=colorsType) +
  theme_light() +
  theme(axis.title = element_text(size=5, face = "bold"),axis.text=element_text(size=3),legend.title=element_text(size=6, face = "bold"),legend.text=element_text(size=4),title=element_text(size=7, face = "bold"),
        legend.spacing= unit(0.1,"points"),legend.spacing.y=unit(0.01,"mm"),legend.key.size = unit(3, "mm"),panel.grid =element_line(color="white"))+
  xlab(paste0("PC1: ", pcavst3.perc.P1,"% variance" )) +
  ylab(paste0("PC2: ", pcavst3.perc.P2,"% variance" )) +
  ggtitle("The first two components of the PCA on whole dataset")
facet1

lgdType_grob=grid.grabExpr(draw(lgdType)) 
facet1_ld<-grid.arrange(facet1,lgdType_grob, ncol=2,widths=c(5,0.8))
facet1_ld

save_plot(paste0(ofig,"PCA_grouping_by_age_time_typecell.png"),facet1_ld)

facets2 <- ggplot(pca_dfA, aes(x=PC1, y=PC2, color=time,fill=time, shape=age)) +
  geom_point(size=0.8,aes(shape=age), color="grey") + 
  scale_fill_manual(values=colorsTime) +
  scale_color_manual(values=colorsTime) +
  scale_shape_manual(values=c(21,25))+
  stat_ellipse(aes(group=time),size=0.2) + 
  theme_light() +
  theme(axis.title = element_text(size=5, face = "bold"),axis.text=element_text(size=3),legend.title=element_text(size=6, face = "bold"),legend.text=element_text(size=4),title=element_text(size=6, face = "bold"),
        legend.spacing= unit(0.1,"points"),legend.spacing.y=unit(0.01,"mm"),strip.text.x = element_text(size = 4),strip.text.y = element_text(size=4),legend.key.size = unit(3, "mm"))+
  xlab(paste0("PC1: ", perc.P1[1],"% variance" )) +
  ylab(paste0("PC2: ", perc.P2[2],"% variance" )) +
  facet_grid(vars(type2))

facets2
save_plot(paste0(ofig,"PCA_typecell_grouping_age_time.png"),facets2, base_width=8, base_height= 3)

facets2.5 <- ggplot(pca_dfA, aes(x=PC1, y=PC2, color=time)) +
  geom_point(size=0.8, show.legend = F) + 
  scale_color_manual(values=colorsTime) +
  scale_x_continuous(n.breaks = 3)+
  stat_ellipse(aes(group=time),size=0.2, show.legend = F) + 
  theme_bw() +
  theme(axis.title = element_text(size=5, face = "bold"),axis.text=element_text(size=3),legend.title=element_text(size=6, face = "bold"),legend.text=element_text(size=4),title=element_text(size=7, face = "bold"),
        legend.spacing= unit(0.1,"points"),legend.spacing.y=unit(0.01,"points"),strip.text.x = element_text(size = 4, face = "bold", color="white"),strip.text.y = element_text(size=4, face = "bold", color="white"),legend.key.size = unit(3, "mm"),panel.spacing = unit(0.5,"mm"),panel.grid =element_line(color="white"))+
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
facet2.5_ld<-grid.arrange(g2.5,lgdTime_grob, ncol=2,widths=c(5,0.8))
facet2.5_ld

pdf(paste0(ofig,"pca_bytype.pdf"), height=4, width=12)
facets2.5
dev.off()

# ==== time splitted:
pdf(paste0(ofig,"pca_bytime.pdf"))
facet3<- ggplot(pca_dfB, aes(x=PC1,y=PC2, color=type, shape=age)) +
  geom_point(size=0.8) +
  scale_color_manual(values=colorsType) +
  stat_ellipse(aes(group=type),size=0.2) + 
  theme_light() +
  theme(axis.title = element_text(size=5, face = "bold"),axis.text=element_text(size=3),legend.title=element_text(size=5, face = "bold"),legend.text=element_text(size=4),title=element_text(size=7, face = "bold"),
        legend.spacing= unit(0.1,"points"),legend.spacing.y=unit(0.01,"points"),strip.text.x = element_text(size = 4),strip.text.y = element_text(size=4),legend.key.size = unit(3, "mm"))+
  facet_grid(. ~ time) +
  ggtitle("PCA on vst sub-matrices at each time")
facet3
dev.off()

# ==== all (time and type) splitted
pdf(paste0(ofig,"pca_typebytime.pdf"))
facet4<-ggplot(pca_dfC, aes(x=PC1,y=PC2, color=age)) +
  geom_point(size=0.2, show.legend = F) +
  scale_color_manual(values = colorsAge) +
  facet_grid_blank(vars(type),vars(time), drop = FALSE) +
  theme_bw() +
  theme(axis.title = element_text(size=5, face = "bold"),axis.text=element_text(size=3),legend.title=element_text(size=4, face = "bold"),legend.text=element_text(size=2),title=element_text(size=7, face = "bold"),
        legend.spacing= unit(0.1,"points"),legend.spacing.y=unit(0.01,"points"),strip.text.x = element_text(size = 4, face = "bold", color="white"),strip.text.y = element_text(size=4, face = "bold", color="white"),legend.key.size = unit(3, "mm"),panel.spacing = unit(0.5,"mm"),panel.grid =element_line(color="white"))+ 
  geom_mark_ellipse(expand = unit(0.8, "mm"),size=0.2, show.legend = F)+
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
lgdAge_grob=grid.grabExpr(draw(left=lgdAge)) 
facet4_ld<-grid.arrange(g4,lgdAge_grob, ncol=2,widths=c(5,0.8))
facet4_ld



save_plot(paste0(ofig,"pca_typebytime.png"),facet4, base_width=10, base_height= 6)

dev.off()
png(paste0(ofig,"PCA_AllData_byCelltype_byCellTypeTime.png"),units = "in", width=10, height= 6, res = 300, family = "Arial")
facet1 / ( facets2.5 |facet4 )
dev.off()
png(paste0(ofig,"PCA_AllData_byCelltype_byCellTypeTime2.png"),units = "in", width=4, height= 5.5, res = 300, family = "Arial")
pcas_plot<-ggpubr::ggarrange(facet1_ld , facet2.5_ld , facet4_ld, ncol=1,heights = c(0.25,0.25,0.5), labels = c("B","C","D") ) 
pca_lgd<-ggpubr::ggarrange(lgdType_grob , lgdTime_grob , lgdAge_grob, ncol=1, heights = c(0.25,0.25,0.5)) 
pcas_plot<-ggpubr::ggarrange(facet1 , lgdType_grob , g2.5, lgdTime_grob , g4, lgdAge_grob, ncol=2, heights = c(0.25,0.25,0.5), widths=c(5,0.8)) 
dev.off()


png(paste0(ofig,"Figure1_head.png"),units = "in", width=8, height= 5.5, res = 300, family = "Arial")
CovHeatmap_grob=grid.grabExpr(draw(heatmap_cor,heatmap_legend_side="bottom", annotation_legend_list=pd2, annotation_legend_side="right",
                                   legend_grouping = "original")) 
ggpubr::ggarrange(CovHeatmap_grob,ggpubr::ggarrange(facet1_ld , facet2.5_ld , facet4_ld, ncol=1,heights = c(0.25,0.25,0.5), labels = c("B","C","D") ) , ncol=2, labels = "A")
dev.off()


#  ====================== DISPERSIONS, SIZE FACTORS =========================
plotDispersionsTissues <- function( fmat, metadata, tisscol = "type",
                                    outfilename, nbco, nbro){
  tissues = sort(unique(metadata[[tisscol]]))
  pdf(outfilename, paper="a4")
  par(mfrow=c(nbro,nbco))
  for (ct in tissues){
    mat.ct <- fmat[,str_detect(colnames(fmat), ct)]
    meta.ct <- metadata %>% filter(str_detect(rownames(metadata),ct))
    x.keep <- apply(mat.ct, 1, function(row) ifelse(sum(row >=5)>= 3, TRUE, FALSE))
    mat.ct <- mat.ct[x.keep,]
    ds.o <- DESeqDataSetFromMatrix(mat.ct, meta.ct, design = ~ age)
    ds.o$age <- relevel(ds.o$age, ref="Young")
    d <- DESeq2::estimateSizeFactors(ds.o)
    d <- DESeq2::estimateDispersions(d)
    plotDispEsts(d, cex=0.6, genecol=rgb(0.1,0.1,0.1,0.3), fitcol="gold", 
                 finalcol=rgb(0.1,0.8,0.7,0.6), main=ct)
  } # end for
  dev.off()
  par(mfrow=c(1,1))
}
 plotDispersionsTissues(fmat, metadata, "type", 
                        outfilename=paste0(ofig,"dispersionsPlot.pdf"), 2, 3)

# ==
# END 
