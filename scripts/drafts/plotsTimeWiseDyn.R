
# import packages
# ============================================================ 
library(dplyr)
library(tidyverse)
library(ggplot2)
library(cowplot)
library(forcats)
library(RColorBrewer)
library(MASS)
library(pheatmap)
library(ggsci) # publishing palettes
library(cowplot)
library(gridExtra)
library(reshape2)
library(ggthemes)
library(scales) # 'viridis_pal' et autres pal
library(ggrepel) # for labels to points
library(org.Mm.eg.db)  
library(AnnotationDbi)
library(edgeR)
# ============================================================ 
# 
setwd("~/BulkAnalysis_plusNetwork")
resdir = "test/"
vrs="vB"
prefil_cou <- "data/prefiltered_counts.rds"
metadata.rds <- "data/metadata.rds"

ens2entrez <- function(vectorensemblids){
  res <- AnnotationDbi::select(org.Mm.eg.db,
                               key=vectorensemblids,
                               columns="ENTREZID", 
                               keytype = "ENSEMBL")
  return(res) # a dataframe columns ENSEMBL and ENTREZID
}
entrez2symbol <- function(vectorensemblids){
  res <- AnnotationDbi::select(org.Mm.eg.db,
                               key=vectorensemblids,
                               columns="SYMBOL", 
                               keytype = "ENTREZID")
  return(res) # a dataframe columns ENSEMBL and ENTREZID
}

# ========================  doing in form of heat map, only for M1, in Old

allkeggenes = getGeneKEGGLinks(species.KEGG = "mmu", convert = FALSE)
dim(allkeggenes) # 36015     2
gokegg <- readRDS( file=paste0(resdir,"edger_dynINTRAkeggo","Old",vrs,".rds"))
DYNgenes <- read.table(paste0(resdir,"edger_dynINTRA_","Old",vrs,".txt"), sep='\t',
                       header=T)
DYNkegg_df <- allkeggenes[allkeggenes$GeneID %in% unique(DYNgenes$id),]
dim(DYNkegg_df) # [1] 2032 2  (only including our more dynamic genes )


names(gokegg[["M1_kegg"]])
trast <- "D4vsD2"
tmpkegg_t_t <- gokegg[["M1_kegg"]][[trast]] %>% 
  filter(N >= mean(gokegg[["M1_kegg"]][[trast]]$N) &
         (P.Up <= 0.0005 | P.Down <=0.0005 )  )  #  has unique pathway ids as rownames ('path:mmu04933')
# addd annotation to add to plot:
tmpkegg_t_t <- tmpkegg_t_t %>% mutate(path_trend = case_when(
  P.Up == P.Down ~ "unknown",
  P.Up > P.Down ~ "Up",
  P.Down > P.Up ~ "Down"
))

localgenes <- DYNgenes %>% filter(type=="M1" & contrast==trast) %>%
  arrange(logFC) %>% mutate(direction = ifelse(logFC > 0, "Up", "Down"))
DYNkegg_genes <- allkeggenes %>% filter(PathwayID %in% rownames(tmpkegg_t_t) &
                                      GeneID %in% unique(localgenes$id))

matrixlfc <- array(0,dim=c(length(unique(DYNkegg_genes$PathwayID)),
  length(unique(localgenes$id))  ) )
rownames(matrixlfc) <- unique(DYNkegg_genes$PathwayID)
colnames(matrixlfc) <- unique(localgenes$id)
for (gene in localgenes$id){
  #print(localgenes[localgenes$id==gene,]$logFC)
  thesepaths <- DYNkegg_genes[DYNkegg_genes$GeneID == gene,] 
  for (path in thesepaths$PathwayID){
    matrixlfc[path, as.character(gene)] <- localgenes[localgenes$id==gene,]$logFC
  }
}


# transform column and rows to names useful to eyes:
labeledm <- matrixlfc
colnames(labeledm) <- entrez2symbol(colnames(matrixlfc))$SYMBOL
rownames(labeledm) <- tmpkegg_t_t[rownames(matrixlfc),]$Pathway
# pulling  the annotations:
#  
direcgenes <- localgenes$direction
names(direcgenes) <- colnames(labeledm)
annotpath <- tmpkegg_t_t %>% dplyr::select(path_trend)
annotpath <- data.frame("path_trend" = annotpath[rownames(matrixlfc),])
rownames(annotpath) <- rownames(labeledm)
# remove zero columns
keep <- apply(labeledm,2, function(x) sum(x)==0)
names(keep) <- as.character(names(keep))
labeledm <-  labeledm[,!keep]

# aesthetics
paletteLength <- 100
myColor <- colorRampPalette(rev(brewer.pal(n = 7, name = "PuOr")))(paletteLength)
# length(breaks) == length(paletteLength) + 1
# use floor and ceiling to deal with even/odd length pallettelengths
myBreaks <- c(seq(min(localgenes$logFC), 0, length.out=ceiling(paletteLength/2) + 1), 
              seq(max(localgenes$logFC)/paletteLength, max(localgenes$logFC), 
                  length.out=floor(paletteLength/2)))



pdf(paste0(resdir,"dynamicTOPgenes_M1_Old_vB.pdf"), width=14, height=7)
pheatmap(labeledm, border_color="white",
         cluster_cols=F, cluster_rows = F,
         color =myColor,
         breaks=myBreaks,
         fontsize = 7,
         angle_col = 45,
         annotation_row = annotpath,
         annotation_colors = list("path_trend" = c(
           "Up" = "indianred", "Down" = "cadetblue", "unknown"="gray")),
         main = "M1 : D4vsD2") + theme_bw()
dev.off()

#### ===== if time add cool color to text
#colorgenes <- rep(c("brown","purple"), length(rownames(labeledm))+1)
p <- pheatmap(labeledm, border_color="white",
              cluster_cols=F, cluster_rows = F,
              color =myColor,
              breaks=myBreaks,
              fontsize = 7,
              angle_col = 45,
              annotation_row = annotpath,
              annotation_colors = list("path_trend" = c(
                "Up" = "indianred", "Down" = "cadetblue", "unknown"="gray")),
              main = "M1 : D4vsD2") 
textcolors <- sapply(annotpath$path_trend, function(x)
  ifelse(x=="Up", "indianred",ifelse(x=="Down", "cadetblue","gray")))
my_gtable = p$gtable
my_gtable$grobs[[4]]$gp=grid::gpar(col=textcolors, fontsize=9)
pdf(paste0(resdir,"dynamicTOPgenes_M1_Old_vC.pdf"),width=15, height=7)
plot(my_gtable) + theme_bw()
dev.off()

# ========================  doing in form of Alluvial

library(ggalluvial) ## include in pkg install 
library(purrr)  ##include in pkg install 

options(stringsAsFactors = F)

#Generate dummy data
myData = map_df(paste0("Person", 1:9), function(x){
  data.frame(     Person = x, 
                  Characteristic = paste0("Characteristic ",
                                          sample(1:5, sample(2:5, 1))))
})
myData$other <- rep(c("otherC","otherB","otherA"), (dim(myData)[1]/3) )
myData$lfc <- rep(c(-2,-1,0,2,3,1.5, 2.2), ((dim(myData)[1]/6)))[1:33]
#https://cran.r-project.org/web/packages/ggalluvial/vignettes/labels.html
#Create the alluvial plot
ggplot(data = myData,
       aes(axis1 = other, axis2=Person, axis3 = Characteristic, y = 0.5)) +
  scale_x_discrete(expand = c(.4, 0)) +
  geom_alluvium(aes(fill = lfc), show.legend = FALSE) +
  geom_stratum() + geom_text(aes(label = Person),
                             stat = "stratum", size = 3) +
  theme_void() + ggtitle("Alluvial plot") + 
  scale_fill_distiller(type = 'seq', palette = 'GnBu',
                       direction = -1, 
                       labels = comma, 
                       guide = "colourbar", 
                       aesthetics = "fill") + 
  ggrepel::geom_text_repel(
    aes(label =  as.character(other)),
    stat = "stratum", size = 4, direction = "y", nudge_x = -.5
  ) + 
  ggrepel::geom_text_repel(
    aes(label =  as.character(Characteristic)),
    stat = "stratum", size = 4, direction = "y", nudge_x = .5
  )
#guides(fill = guide_colorbar(barwidth = 0.5, barheight = 10, title="Value"))


ggplot(myData,
       aes(x = other, stratum = Person, alluvium = Characteristic, y = 0.5,
           fill = response)) +
  scale_x_discrete(expand = c(.4, 0)) +
  geom_flow(width = 1/4) +
  geom_stratum(alpha = .5, width = 1/4) +
  scale_linetype_manual(values = c("blank", "solid")) +
  ggrepel::geom_text_repel(
    aes(label = ifelse(as.numeric(survey) == 1, as.character(response), NA)),
    stat = "stratum", size = 4, direction = "y", nudge_x = -.5
  ) +
  ggrepel::geom_text_repel(
    aes(label = ifelse(as.numeric(survey) == 3, as.character(response), NA)),
    stat = "stratum", size = 4, direction = "y", nudge_x = .5
  ) +
  theme(legend.position = "none") +
  ggtitle("vaccination survey responses", "labeled using `geom_text_repel()`")

