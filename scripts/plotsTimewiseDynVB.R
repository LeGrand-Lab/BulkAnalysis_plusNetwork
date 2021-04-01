
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

# set locations
# ============================================================ 

setwd("~/BulkAnalysis_plusNetwork")
resdir = "plotsDE/"
DEdir = "signaturestypes/"
prefil_cou <- "data/prefiltered_counts.rds"
metadata.rds <- "data/metadata.rds"
# ============================================================ 
#functions
# ============================================================ 
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
# ============================================================ 

#set vars
vrs="vB"
age= "Young"
# ============================================================ 
# ========================  doing in form of heat map, 

# allkeggenes = getGeneKEGGLinks(species.KEGG = "mmu", convert = FALSE)
# dim(allkeggenes) # 36015     2
gokegg <- readRDS( file=paste0(DEdir,"edger_dynINTRAkeggo",age,vrs,".rds"))
DYNgenes <- read.table(paste0(DEdir,"edger_dynINTRA_",age,vrs,".txt"), sep='\t',
                       header=T)
DYNkegg_df <- allkeggenes[allkeggenes$GeneID %in% unique(DYNgenes$id),]
dim(DYNkegg_df) # [1] 2032 2  (only including our more dynamic genes )

#detect cell types 
cts <- unique(DYNgenes$type) 
names(gokegg[[ki]])

for (ct in cts){
  skip_to_next <- FALSE
    tryCatch({
      print(ct)
      ki <- paste0(ct,"_kegg")
      trastes <- names(gokegg[[ki]])
      for(trast in trastes){
        ## **
        print(dim(gokegg[[ki]][[trast]]))
        # -------------------------- innermost loop 
        # -------------------------------------------------------
        tmpkegg_t_t <- gokegg[[ki]][[trast]] %>% 
          filter(N >= mean(gokegg[[ki]][[trast]]$N) &
                   (P.Up <= 0.0005 | P.Down <=0.0005 )  )  #  has unique pathway ids as rownames ('path:mmu04933')
        # addd annotation to add to plot:
        tmpkegg_t_t <- tmpkegg_t_t %>% mutate(path_trend = case_when(
          P.Up == P.Down ~ "unknown",
          P.Up > P.Down ~ "Up",
          P.Down > P.Up ~ "Down"
        ))
        # -------------------------------------------------------
        # -------------------------- innermost loop 
        # -------------------------------------------------------
        localgenes <- DYNgenes %>% filter(type==ct & contrast==trast) %>%
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
        # -------------------------------------------------------
        # -------------------------- innermost loop 
        # -------------------------------------------------------
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
        # -------------------------------------------------------
        # -------------------------- innermost loop 
        # -------------------------------------------------------
        # aesthetics
        paletteLength <- 100
        myColor <- colorRampPalette(rev(brewer.pal(n = 7, name = "PuOr")))(paletteLength)
        # length(breaks) == length(paletteLength) + 1
        # use floor and ceiling to deal with even/odd length pallettelengths
        myBreaks <- c(seq(min(localgenes$logFC), 0, length.out=ceiling(paletteLength/2) + 1), 
                      seq(max(localgenes$logFC)/paletteLength, max(localgenes$logFC), 
                          length.out=floor(paletteLength/2)))
        theme_set(theme_bw())
        p <- pheatmap(labeledm, border_color="white",
                      cluster_cols=F, cluster_rows = F,
                      color =myColor,
                      breaks=myBreaks,
                      fontsize = 7,
                      angle_col = 45,
                      annotation_row = annotpath,
                      annotation_colors = list("path_trend" = c(
                        "Up" = "indianred", "Down" = "cadetblue", "unknown"="gray")),
                      main = paste0(age," ",ct," ",trast )) 
        textcolors <- sapply(annotpath$path_trend, function(x)
          ifelse(x=="Up", "indianred",ifelse(x=="Down", "cadetblue","gray")))
        my_gtable = p$gtable
        my_gtable$grobs[[4]]$gp=grid::gpar(col=textcolors, fontsize=9)
        # -------------------------------------------------------
        # -------------------------- innermost loop 
        # -------------------------------------------------------
        if (ct %in% c("FAPs","ECs","sCs")){
          pdf(paste0(resdir,"dynamicTOPgenes_",age, ".", ct, ".",trast, ".",vrs,".pdf"),width=10, height=4)
          plot(my_gtable) 
          dev.off()
        } else{
          pdf(paste0(resdir,"dynamicTOPgenes_",age, ".", ct, ".",trast, ".",vrs,".pdf"),width=15, height=7)
          plot(my_gtable) 
          dev.off()}
        
        # -------------------------------------------------------
        # ---- end innermost loop
      }
        ## **
  },error=function(e){skip_to_next <<- TRUE})#end TryCatch
  if(skip_to_next==T){next}
}
  