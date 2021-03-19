library(dplyr)
library(tidyverse)
library(ggplot2)
library(cowplot)
library(forcats)
library(RColorBrewer)
library(MASS)
library(pheatmap)
library(DESeq2)
library(ggsci) # publishing palettes
library(cowplot)
library(gridExtra)
# devtools::install_github("vqv/ggbipplot")
# install.packages("glmpca")
setwd("~/bulk_analysis/")

ofig <- "plotsPrelim/"
#outputs:
prefil_cou <- "data/prefiltered_counts.rds"
prefil_tpm <- "data/prefiltered_TPM.rds"
metadata.rds <- "data/metadata.rds"

# REDO here: enlever AllNEg, et aussi les REDO
# preliminary counts and TPM matrices formatting
# =========================================================================
pc.mat <- readRDS("data/protcod_counts.rds")
pc.TPM <- readRDS("data/protcod_TPM.rds")
design <- read.table("data/design.csv",header=T,sep=",", row.names=1)

# filter out all-zeroes count rows:
fTPM <- pc.TPM[!rowSums(pc.mat)==0,]  #c'est bien fmat criterium, MUST first!
fmat <- pc.mat[!rowSums(pc.mat)==0,]  

fTPM <- as.matrix(fTPM, dim=dim(fTPM))
fmat <- as.matrix(fmat, dim=dim(fmat))
typeof(fmat[3000,50])
#[1] "integer"
absentsamp <- setdiff(colnames(fmat),design$sample)
# samples not real ones:
#[1] "X2.422363" "X2.422364" "X2.424030" "X2.424031" "X2.428151" "X2.428100" "X2.429794" "X2.429793"
exclu <- fmat[,absentsamp]
fmat <- fmat[,!(colnames(fmat) %in% absentsamp)]
fTPM <- fTPM[,!(colnames(fTPM) %in% absentsamp)]

# take away REDO and All_Neg from matrices and metadata
metadata <- design
metadata <- metadata[!str_detect(metadata$type,"All_Neg"),]
metadata <- metadata[!str_detect(metadata$type,"REDO"),]
#retain only sample names present in metadata:
fmat <- fmat[,colnames(fmat) %in% metadata$sample]
fTPM <- fTPM[,colnames(fTPM) %in% metadata$sample]

# rename sample names:
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
# dfnwn <- dfnwn %>% mutate(
#     newnames = str_replace(
#       str_replace(
#         paste0(str_replace(nwn,"-RNA",""),"_",replicate),
#       "Old","O"),
#       "Young","Y"))
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

# save these new matrices, and metadata
saveRDS(fTPM, file=prefil_tpm)
saveRDS(fmat, file=prefil_cou)
saveRDS(metadata, file=metadata.rds)

# =========================================================================
fmat <- readRDS(prefil_cou)
fTPM <- readRDS(prefil_tpm)
metadata <- readRDS(metadata.rds)
##
?cor
cor.mat <- cor(fTPM,method="spearman")

pheatmap(cor.mat, fontsize = 7)

# TODO: add barcolors to visualize celltypes and age
#  testing same cor function on log2(TPM+1) data : same result
cormatB <- cor(log(fTPM+1))
pheatmap(cormatB, method="pearson",fontsize = 7)
# =========================================================================
# make DESeq2 object
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
  
pcavst2 <- plotPCA(vst, intgroup=c("time","age"),returnData=T)
gg2 <- ggplot(pcavst2, aes(x=PC1, y=PC2, color=time, shape=age)) +
  geom_point(size=1) + theme_light() + stat_ellipse() +
  scale_color_manual(values=mycol) +
  xlab(paste0("PC1: ", percentVar[1],"% variance" )) +
  ylab(paste0("PC2: ", percentVar[2],"% variance" )) +
  ggtitle("", subtitle="grouping by age + time")

pcavst3 <- plotPCA(vst, intgroup=c("type", "time","age"), returnData=T)
gg3 <- ggplot(pcavst3, aes(x=PC1, y=PC2,  color=type, shape=time)) +
  geom_point(size=2) + 
  stat_ellipse(aes(group=type)) + 
  scale_color_manual(values=rev(mycol)) +
  theme_light() +  
  xlab(paste0("PC1: ", percentVar[1],"% variance" )) +
  ylab(paste0("PC2: ", percentVar[2],"% variance" )) +
  ggtitle("PCA when grouping by type + time + age")

pdf(paste0(ofig, "samplesPCA_2.pdf"))
plot_grid(plot_grid(gg1,gg2,nrow=1),gg3, nrow=2,rel_heights=c(1,2))
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

facets2 <- ggplot(pca_dfA, aes(x=PC1, y=PC2, color=time, shape=age)) +
  geom_point(size=2, alpha=.6) + 
  scale_color_d3()+
  theme_light() +
  xlab(paste0("PC1: ", perc.P1[1],"% variance" )) +
  ylab(paste0("PC2: ", perc.P2[2],"% variance" )) +
  facet_grid(. ~ type2)

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
ggplot(pca_dfC, aes(x=PC1,y=PC2, color=age)) +
  geom_point(size=2, alpha=.6) +
  scale_color_d3() + theme_light() +
  facet_grid(vars(time),vars(type)) +
  ggtitle("PCA on vst sub-matrices time_type")
dev.off()
# ==

# ==  using rlog (just for comparison vs. vst  )
# rlog() may take a long time with 50 or more samples,
#  vst() is a much faster transformation
# rld <- rlog(ds.o) --> NOT done as TIME CONSUMING 
# pcavst.rld <- plotPCA(rld, intgroup=c("time","age","type"),returnData=T)



# ggplot(gpca.dat, aes(x = PC1, y = PC2, color = dex, shape = cell)) +
#  geom_point(size =3) + coord_fixed() 


####################
# pc.vst <- prcomp(assay(vst))
# pc.vst.df <- as.data.frame(pc.vst$rotation[,2:3])
# ggplot(pc.vst.df, aes(x = PC2, y = PC3)) +
#   geom_point(size =1) 

# set.seed(123)
# trmat <- as.matrix(fmat[sample(1:20000,2000),]) # take tiny  amount of genes
# pca <- prcomp(trmat, center=T, scale=T) 
# library(ggbiplot)
# ggbiplot(pca)
# library("glmpca")
# gpca <- glmpca(trmat, L=5)
# gpca.df <- gpca$factors
# gpca.df$celltype <- 
#   gpca.df$timepoint <-
#   gpca.df$age <-
#   gpcax <- glmpca(trmat, L=2)



# matcov <- cov(fmat[sample(1:21000,2000),]) 
# scaledcov <- diag(sqrt(1/diag(matcov))) %*% matcov
# m_eigen <- eigen(matcov)
# mm_eigen <- eigen(scaledcov)

# mycol <- c("gray","lightblue","blue","purple","pink","gold","darkgreen")


####
# pcavst4 <- plotPCA(vst, intgroup=c("group"), returnData=T)
# ggplot(pcavst4, aes(x=PC1, y=PC2, color=group)) +
#   geom_point(size=2) + theme_light() + stat_ellipse()
# # ==== 2 attempt
# ds.o2 <- DESeqDataSetFromMatrix(countData=fmat,
#                                 colData=metadata,
#                                 design= ~type+time)
# 
# vst2 <- varianceStabilizingTransformation(ds.o2)
# 
# pcavst1.2 <- plotPCA(vst2, intgroup=c("time"),returnData=T)
# ggplot(pcavst1.2, aes(x=PC1, y=PC2, color=time)) +
#   geom_point(size=1) + theme_light()
# 
# pcavst2.2 <- plotPCA(vst2, intgroup=c("time","type"),returnData=T)
# ggplot(pcavst2.2, aes(x=PC1, y=PC2, color=type, shape=time)) +
#   geom_point(size=3,alpha=.3) + theme_light() + stat_ellipse()
#######
