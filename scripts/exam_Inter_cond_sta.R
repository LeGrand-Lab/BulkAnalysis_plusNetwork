## 
# Detect Differential Expression between Young and Old at each day-point, by cellType
# output: all results into 'exam_INTER_conditions/static/'
# --
# Joha GL + MoulleP
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
# Working directory
setwd("~/BulkAnalysis_plusNetwork/")
# Raw data directory
inData<- "data/"
# Save count normalided in data directory
NormData <- "data/CountNormalised/"
# Directory of static DEG analysis
odir <- "exam_INTER_conditions/"
staticD = paste0(odir,"static/")

#Created directories
listDirectorie<-c(inData,NormData,odir,staticD)
for ( directories in listDirectorie){
  if (dir.exists(directories) == F) {
    dir.create(directories)
  }
}
listInterCondtionDirectories<-c(staticD)
for ( sd in c("TableDEG/","PlotsDEG/","analysis/") ){
  for (d in listInterCondtionDirectories ){
    if (dir.exists(paste0(d,sd)) == F) {
      dir.create(paste0(d,sd))
    }
  }
}

#Loads 
# Correspondance geneID and symbol
genes_df <- read.table("data/genesinfo.csv", sep="\t", header=T)
#Counting matrix after first filtering
prefil_cou <- paste0(inData,"prefiltered_counts.rds")
#Infos about samples
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
orderTypecell=c("ECs","FAPs","MuSCs","Neutrophils","Inflammatory-Mac","Resolving-Mac")

colorsType=c("#10b387ff","#3d85c6ff","#b171f1ff","#f0e442ff","#ff9900ff","#cc0000ff")
names(colorsType)=orderTypecell
colorsTime = c(brewer.pal(9,"BuPu")[3],brewer.pal(9,"BuPu")[5],brewer.pal(9,"BuPu")[7],brewer.pal(9,"BuPu")[9])


## Load files
fmat <- readRDS(prefil_cou)
metadata <- readRDS(metadata.rds)
metadata <- metadata %>% mutate(day=time) %>% mutate(daytype=paste0(day,".",type)) 
#keep genes with expression > 5 on at least 3 samples
keep <- apply(fmat, 1, function(row) ifelse(length( which( row >= 5)) >= 3, TRUE, FALSE) )
fmat <- fmat[keep,]


#####
# DE analysis
# static old vs young
# by cell-type by day
#####

##  with 'subset' s   on same original object:
#  ===========================================================================
deseqdataset.static <- DESeqDataSetFromMatrix(countData = fmat,
                               colData = metadata,
                               design= ~ age)
# set "Young" as reference
colData(deseqdataset.static)$age <- factor(colData(deseqdataset.static)$age,
                            levels=c("Young", "Old"))
#initialization result DEG analyse to concatenated results
day.type <- sort(unique(colData(deseqdataset.static)$daytype))

sortie_full <- data.frame("baseMean"= double(), "log2FoldChange"= double(), 
                          "lfcSE"= double(),  "stat" = double(), 
                          "pvalue" = double(), "padj" = double(),
                          "id" = character(), "day" = character(), "type"=character())
out_softfilt<-data.frame()


Firstday<-str_split(day.type[1], "\\.")[[1]][1]
# Initialization Total table normalized Young and Old count
  #mean replicat sample 
  countNormalisedtable<-data_frame(geneID = character(),age=character(),celltype=character(),day=character(),meancountnormalised=character())
  #keep replicats
  countNormalised1<-data.frame("gene"=rownames(fmat))
for (k in day.type){
  day <- str_split(k, "\\.")[[1]][1]
  type <- str_split(k, "\\.")[[1]][2]
  print(paste(day,type))
  # do subsets as needed + deseq analysis
    deseqdataset.static.daytype <- deseqdataset.static[,deseqdataset.static$type == type & deseqdataset.static$day == day]
    deseq.static.daytype <- DESeq(deseqdataset.static.daytype, test="Wald", full=~age)
  
  # Get countNormalised 
    countNormalised<-counts(deseq.static.daytype , normalize =T)
    countNormalised1<-cbind(countNormalised1,countNormalised)
    
    #mean replicats
    countNormalisedtable_tempo= data.frame(geneID=c(rownames(countNormalised),rownames(countNormalised)),
                                           age=c(rep("Young", dim(countNormalised)[1]),rep("Old", dim(countNormalised)[1])),
                                           celltype=rep(type, dim(countNormalised)[1]),
                                           day=rep(day, dim(countNormalised)[1]),
                                           meancountnormalised=c(rowMeans(countNormalised[,str_detect(colnames(countNormalised),"Young")]),rowMeans(countNormalised[,str_detect(colnames(countNormalised),"Old")])))
    
    countNormalisedtable<- rbind(countNormalisedtable,countNormalisedtable_tempo)

  # Add Shrink log2 fold changes
    deres <- lfcShrink(deseq.static.daytype, coef="age_Old_vs_Young", type="apeglm",
                       parallel=T)
    deres$id <- rownames(deres)
    deres$day = day
    deres$type = type
  # solft filters
    sortie_full <- rbind(sortie_full, as_tibble(deres))
    
    deres_f <- as_tibble(deres) %>% dplyr::filter(padj <= 0.5 & abs(log2FoldChange) >= 0.1)
    out_softfilt <- rbind(out_softfilt,deres_f)
}
# Finalization Get countNormalised 
countNormalisedtable$symbol<-genes_df[match(countNormalisedtable$geneID,genes_df$Geneid),]$symbol
countNormalised1$symbol<-genes_df[match(rownames(countNormalised1),genes_df$Geneid),]$symbol
write.table(countNormalisedtable,paste0(NormData,"MeanCountNormalised.txt"),row.names = F, sep = '\t')
write.table(countNormalised1,paste0(NormData,"CountNormalised.txt"),row.names = F, sep = '\t')


# complete tables with gene_symbols & save table
sortie_full$symbol <- genes_df[match(sortie_full$id, genes_df$Geneid),]$symbol
sortie_full<- sortie_full %>% dplyr::filter(baseMean!=0) 
saveRDS(sortie_full, paste0(staticD, DEG_table, "static_full.rds")) 
write.table(sortie_full, paste0(staticD,DEG_table, "static_full.csv"), sep="\t", col.names=T, row.names = F) 
out_softfilt <- sortie_full %>% dplyr::filter(padj <= 0.5 & abs(log2FoldChange) >= 0.1)
write.table(out_softfilt, paste0(staticD, DEG_table, "static_softfilter.csv"), sep="\t", col.names=T, row.names = F) 
out_hardfilt = out_softfilt %>% dplyr::filter(padj <= 0.05 & abs(log2FoldChange) >= 1.2)
write.table(out_hardfilt, paste0(staticD, DEG_table, "static_hardfilter.csv"), sep="\t", col.names=T, row.names = F)



# ** VOLCANOES ** :
if (volcanoneeded){
 
  sortie_full <- read.table(paste0(staticD, DEG_table, "static_full.csv") ,sep="\t", header=T)
  names(colorsTime) = sort(unique(sortie_full$day))
  sortie_full <- sortie_full %>% mutate (type = factor( sortie_full$type , levels =orderTypecell ))
  
 # set aesthetics data:
    sortie_full <- sortie_full %>% mutate(DEsignificant=case_when(
    padj <= 0.05 & log2FoldChange >= 1.2 ~ "Signi UP" ,
    padj <= 0.05 & log2FoldChange <= -1.2 ~ "Signi DOWN",
    TRUE  ~ "Not Signi"
    ))

  g <-ggplot(sortie_full, aes(x=log2FoldChange, y = -log10(padj+1e-20),color=DEsignificant)) +
    geom_point(aes(color=DEsignificant),size=.3) +
    scale_color_manual(values=c("lightgrey",brewer.pal(10, "RdBu")[9],brewer.pal(10, "RdBu")[2] )) +
    geom_vline(xintercept = c(1.2,-1.2), data=,color= "black",
               linetype="dashed", size=.2) +
    facet_grid_blank(vars(day),vars(type), drop = FALSE) + 
    theme_calc() +
    theme(panel.grid.major=element_blank()) +
    
    geom_text_repel(
      data= subset(sortie_full, padj < 0.0005 & 
                     abs(log2FoldChange) > 1.2),
      aes(label=symbol, fill=DEsignificant),
      size=2,
      segment.size = .1,
      force=2, force_pull = 2,
      max.overlaps=15
    ) +
    labs(title="Old vs Young across day&type", subtitle = "Signi = Differential gene padj < 0.05, -1.2<log2FoldChange>1.2",
         caption="vertical lines: ABS(log2FoldChange)=1.2 
          labels only for genes padj < 0.0005")+
    theme(axis.title = element_text(size=5, face = "bold"),axis.text=element_text(size=3),legend.title=element_text(size=6, face = "bold"),legend.text=element_text(size=4),title=element_text(size=7, face = "bold"),
          legend.spacing= unit(0.1,"points"),legend.spacing.y=unit(0.01,"points"),strip.text.x = element_text(size = 6, face = "bold", color="white"),strip.text.y = element_text(size=6, face = "bold", color="white"),legend.key.size = unit(3, "mm"),panel.spacing = unit(0.5,"mm"),panel.grid =element_line(color="white"))
  
  #Add color in strip background
  g2 <- ggplot_gtable(ggplot_build(g))
  stripRowName <- which(grepl('strip-', g2$layout$name))
  k <- 1
  fills <- c(colorsType,colorsTime)
  for (i in stripRowName) {
    j <- which(grepl('rect', g2$grobs[[i]]$grobs[[1]]$childrenOrder))
    g2$grobs[[i]]$grobs[[1]]$children[[j]]$gp$fill <- fills[k]
    k <- k+1
  }
  png(paste0(staticD,"PlotsDEG/volcano_static.png"),units = "in", width=10, height= 5.5, res = 300, family = "Arial")
  grid::grid.draw(g2)
  dev.off()
  tiff(paste0(staticD,"PlotsDEG/volcano_static.tiff"),units = "in", width=10, height= 4, res = 300, family = "Arial")
  grid::grid.draw(g2)
  dev.off()
  pdf(paste0(staticD,"PlotsDEG/volcano_static.pdf"), width=14, height=10)
  grid::grid.draw(g2)
  dev.off()
  
}


