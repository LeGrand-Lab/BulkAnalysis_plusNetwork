##
# I. TPM across libraries Pheatmap. II. RAW COUNTS --> vst --> PCA.
# johaGL 2021
##
library(dplyr)
library(tidyverse)
library(ggplot2)
library(cowplot)
library(RColorBrewer)
library(MASS)
library(pheatmap)
library(DESeq2)
library(ggsci) # publishing palettes
library(cowplot)
library(gridExtra)
library(ggforce)
library("factoextra")
setwd("~/BulkAnalysis_plusNetwork2/")

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
pc.TPM.tempo <- data.frame(pc.TPM ,genes_df[match(rownames(pc.TPM), genes_df$Geneid),]$symbol)
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

batch = read.table("data/batchesinfo.csv") # to add in visuals
metadata$batch = batch[match(metadata$sample, batch$sample),]$batch

myannot = data.frame(CellType=metadata$type,
                           Time = metadata$time,
                     Age = metadata$age,
                     Batch = metadata$batch)
rownames(myannot) = rownames(metadata)

colorsType = brewer.pal(6, "Accent")
names(colorsType) = unique(myannot$CellType)
colorsTime = brewer.pal(4, "Purples")
names(colorsTime) = sort(unique(myannot$Time))
colorsAge = c("black", "gray")
names(colorsAge) = c("Old","Young")
colorsBatch = brewer.pal(4, "Paired")
names(colorsBatch) = sort(unique(myannot$Batch))
mycolors = list("CellType"=colorsType,
              "Time" = colorsTime,
              "Age" =colorsAge , "Batch" = colorsBatch)
pdf(paste0(ofig,"Heatmap_libraries.pdf"), width=10, height= 8 )
heatmap_lib<-pheatmap(cor.mat, fontsize = 8,
    color = colorRampPalette(c("darkgray","white" ,"orange1"))(100),
     cluster_rows = TRUE,
     cluster_cols = TRUE, clustering_distance_rows = "euclidean",
     clustering_distance_cols = "euclidean",  annotation_col = myannot,
     annotation_colors=mycolors, fontsize_row = 5, fontsize_col = 5,
     angle_col = "315", 
    main= "TPM correlation among libraries (Spearman rho scores)")
heatmap_lib
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
ds.o <- DESeqDataSetFromMatrix(countData=fmat,
                                colData=metadata,
                                design= ~type + time + age:time)

# use vst transformation
dim(counts(ds.o))
vst <- varianceStabilizingTransformation(ds.o)
mycol <- pal_jama(alpha=.7)(6)


pcavst1 <- plotPCA(vst, intgroup=c("age"),returnData=T)
percentVar <- round(100* attr(pcavst1,"percentVar"))
gg1 <- ggplot(pcavst1, aes(x=PC1, y=PC2, color=age)) +
  geom_point(size=1) + theme_light() + scale_color_manual(values=mycol) +
  xlab(paste0("PC1: ", percentVar[1],"% variance" )) +
  ylab(paste0("PC2: ", percentVar[2],"% variance" )) +
  ggtitle("PCA", subtitle="grouping by age")
gg1  
save_plot(paste0(ofig,"PCA_grouping_by_age.png"),gg1)

pcavst2 <- plotPCA(vst, intgroup=c("time","age"),returnData=T)
gg2 <- ggplot(pcavst2, aes(x=PC1, y=PC2, color=time, shape=age)) +
  geom_point(size=1) + theme_light() + stat_ellipse() +
  scale_color_manual(values=mycol) +
  xlab(paste0("PC1: ", percentVar[1],"% variance" )) +
  ylab(paste0("PC2: ", percentVar[2],"% variance" )) +
  ggtitle("", subtitle="grouping by age + time")
gg2
save_plot(paste0(ofig,"PCA_grouping_by_age_time.png"),gg2)

pcavst3 <- plotPCA(vst, intgroup=c("type", "time","age"), returnData=T)

gg3 <- ggplot(pcavst3, aes(x=PC1, y=PC2,  color=type, shape=time)) +
  geom_point(size=2) + 
  stat_ellipse(aes(group=type)) + 
  scale_color_manual(values=rev(mycol)) +
  theme_light() +
  xlab(paste0("PC1: ", percentVar[1],"% variance" )) +
  ylab(paste0("PC2: ", percentVar[2],"% variance" )) +
  ggtitle("PCA when grouping by type + time + age")
gg3
save_plot(paste0(ofig,"PCA_grouping_by_age_time_typecell.png"),gg3)
pdf(paste0(ofig, "samplesPCA_2.pdf"))
PCA_grouping<-plot_grid(plot_grid(gg1,gg2,nrow=1),gg3, nrow=2,rel_heights=c(1,2))
PCA_grouping
save_plot(paste0(ofig,"PCA_gridplot_grouping_by_age_time_typecell.png"),PCA_grouping, base_width=12, base_height= 8)

dev.off()

# ====
typescell <- c("FAPs","ECs", "Neutro", "M1", "M2", "sCs")
#metadataA  <- metadata %>% mutate()

pca_dfA <- data.frame("PC1"=numeric(),   
                     "PC2"=numeric(),
                     "group"=character(), "time"=character(),
                     "age" =character(),"name"=character() ,
                     "type2"=character())
perc.P1 <- c()
perc.P2 <- c() 
for (k in typescell){
   subset <- ds.o[,ds.o$type==k]
   subset.vst <- varianceStabilizingTransformation(subset)
   tmp <- plotPCA(subset.vst, intgroup=c("time","age"), returnData=T)
   tmp$type2 <- k
   pca_dfA <- rbind(pca_dfA,tmp)
   percentVar.subset <- round(100* attr(tmp,"percentVar"))
   perc.P1 <- c(perc.P1, percentVar.subset[1])
   perc.P2 <- c(perc.P2, percentVar.subset[2])
}

facets1 <- ggplot(pca_dfA, aes(x=PC1, y=PC2, color=age, shape=time)) +
  geom_point(size=2, alpha=.6) + 
  scale_color_d3()+
  theme_light() +
  xlab(paste0("PC1: ", perc.P1[1],"% variance" )) +
  ylab(paste0("PC2: ", perc.P2[2],"% variance" )) +
  facet_grid(. ~ type2)

facets1
facets2 <- ggplot(pca_dfA, aes(x=PC1, y=PC2, color=time, shape=age)) +
  geom_point(size=2, alpha=.6) + 
  scale_color_d3()+
  theme_light() +
  xlab(paste0("PC1: ", perc.P1[1],"% variance" )) +
  ylab(paste0("PC2: ", perc.P2[2],"% variance" )) +
  facet_grid(. ~ type2)
facets2
save_plot(paste0(ofig,"PCA_typecell_grouping_age_time.png"),facets2, base_width=8, base_height= 3)
plot_grid(facets1, facets2, nrow=2)

facets3 <- ggplot(pca_dfA, aes(x=PC1, y=PC2, color=time, shape=age)) +
  geom_point(size=2, alpha=.6) + 
  scale_color_d3()+
  theme_light() +
  xlab(paste0("PC1: ", perc.P1[1],"% variance" )) +
  ylab(paste0("PC2: ", perc.P2[2],"% variance" )) +
  facet_grid(. ~ type2)

pdf(paste0(ofig,"pca_bytype.pdf"), height=4, width=12)
facets3
dev.off()

# ==== time splitted:
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
  tmp <- plotPCA(subset.vst, intgroup=c("age"), returnData=T)
  pca_dfB <- rbind(pca_dfB,tmp)
  percentVar.subset <- round(100* attr(tmp,"percentVar"))
  pP1 <- c(pP1, percentVar.subset[1])
  pP2 <- c(pP2, percentVar.subset[2])
}
pca_dfB$type <- metadata[match(rownames(metadata),rownames(pca_dfB)),]$type 
pca_dfB$time <- metadata[match(rownames(metadata),rownames(pca_dfB)),]$time 

pdf(paste0(ofig,"pca_bytime.pdf"))
ggplot(pca_dfB, aes(x=PC1,y=PC2, color=age, shape=type)) +
  geom_point(size=2, alpha=.6) +
  scale_color_d3() + theme_light() +
  facet_grid(vars(time)) +
  ggtitle("PCA on vst sub-matrices 'time'")
dev.off()
# ==== all (time and type) splitted
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

pdf(paste0(ofig,"pca_typebytime.pdf"))
pcatypebytime<-ggplot(pca_dfC, aes(x=PC1,y=PC2, color=age)) +
  geom_point(size=2, alpha=.6) +
  scale_color_d3() + theme_light() + geom_mark_ellipse(expand = unit(1, "mm"),)+
  facet_grid(vars(time),vars(type)) +
  ggtitle("PCA on vst sub-matrices time_type")
pcatypebytime
save_plot(paste0(ofig,"pca_typebytime.png"),pcatypebytime, base_width=10, base_height= 6)

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
# plotDispersionsTissues(fmat, metadata, "type", 
#                        outfilename=paste0(ofig,"dispersionsPlot.pdf"), 2, 3)

# ==
# END 
