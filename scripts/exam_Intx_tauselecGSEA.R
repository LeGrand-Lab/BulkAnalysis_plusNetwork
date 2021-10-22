# Test Differential enrichment of pathways
# BY DAY, on unique symbols (unique symbol:celltype) based on Tau filtering
# outputs into 'exam_inter_conditions/static' :
#           rds/Intx_shot_onTauExtract.rds
#           
# JohaGL
setwd("~/BulkAnalysis_plusNetwork/")
library(tidyverse)

library(DESeq2)
library(ComplexHeatmap)
library(circlize)
library(viridis)
library(fgsea)
library(msigdbr)
library(gridExtra)
library(cowplot)
library(RColorBrewer)

M = 1000 # M genes up + M genes down, total of 2M genes for GSEA (each celltype)

odir = "exam_INTER_conditions/static/"
taudir = "Tau/"
consensusfile <- "conseTau_ensemblid.rds"

doDEtest <- F # done
doGSEAintx <- F # done
dopreplots <- F # done
daysv = c("D0","D2","D4","D7")
CUTOFFTAU = 0.3
fmat <- readRDS("data/prefiltered_counts.rds")
metadata <- readRDS("data/metadata.rds")
genes_df <- read.table("data/genesinfo.csv", sep="\t", header=T)

gdir <- "GSEA_Intx/"

print(paste0("cd ",odir, "; if [ ! -d ", 
              gdir,"   ] ;then mkdir -p ", gdir,"rds/; cd ", gdir, "; mkdir csv/; fi"))

# ======== Import tauconsensus and give uniqueness: function importdedup() =======
importdedup <- function(consensusfile, CUTOFFTAU){
    consensus_tau <- readRDS(paste0(taudir, consensusfile))
    # extract only Tau > CUTOFFTAU and deduplicate with max tau
    dedup <- list()
    for (d in daysv){
      tf <- consensus_tau[[d]] %>% filter(Tau >= CUTOFFTAU) %>% 
        group_by(symbol) %>% slice_max(Tau, n=1, with_ties = F)  ### ??? is this ok ????
      dedup[[d]] <- tf
      } # end for
    return(dedup)
}
cat("read how 'consensus' is determined and why, in folder 'exam_INTER_conditions/static/readmedocs' :",
    "\n","file recap_explain.md, section **Tau classification and consensus** ")
# !!  use here function importdedup
dedup <- importdedup(consensusfile, CUTOFFTAU)
for (i in names(dedup)){
  write.table(dedup[[i]] %>% select(symbol, whichMAX ) ,
              paste0(odir, "GSEA_Intx/csv/UNIQUEsym_ct_",i,".csv"), sep='\t', col.names=T)
}

# =============================================================================

if (!doDEtest){
  print("you set 'doDEtest' variable as FALSE, which means already in static/rds the DEG list exists")
}else{print("doDEtest is TRUE, will run DESeq2")}

# ===================== Preliminary heatmaps : vst values day by day ====================================

preliminaryheatmapvst <- function(d, destinyfilepdf ){
  # Z variable is a dataframe of geneid, symbol.x, Tau, whichMAX, nbMAX, exclusiveOld, symbol.y, (where .x is equal to .y))
  Z <- left_join(dedup[[d]], genes_df, by=c("id"="Geneid")) # ensemblids
  Z$symbol <- Z$symbol.x    # one symbol can correspond to two distinct ids 
  Z$symbol.y <- Z$symbol.x <- NULL
  smx <- fmat[rownames(fmat) %in% Z$id, str_detect(colnames(fmat),d)]
  smm <- metadata %>% filter(time==d)
  #  rownames for this day must be unique 
  if (length(unique(rownames(smx))) == length(rownames(smx))){
    print("ok, rownames are unique (ensemblid)")
  }else{
    print("stop, rownames must be set unique")
  }
  
  predsv <- DESeqDataSetFromMatrix(countData=smx, colData=smm, design = ~age)
  dsv <- DESeq2::varianceStabilizingTransformation(predsv)  # dsv <- DESeq2::rlog(predsv)
  mysplitcol = data.frame(x=sapply(colnames(assay(dsv)),
                                   function(x) unlist(str_split(x, "\\."))[2]))
  rownames(mysplitcol) = colnames(assay(dsv))
  # reorder matrix, begin with a empty, then give rownames :
  
  orderrows = Z[match(rownames(assay(dsv)), Z$id),]
  
  givemereorderedMX <- function(amatrix, orderrw, ordercl){
    empty <- array(NA, dim=dim(amatrix))
    colnames(empty) <- ordercl
    rownames(empty) <- orderrw
    for (r in rownames(empty)){
      for (c in colnames(empty)){
        empty[r,c] <- amatrix[r,c]
      }
    }
    return(empty)
  }
  tempocol = data.frame("c" = rownames(mysplitcol),"x"=mysplitcol$x)
  a = orderrows %>% arrange(str_to_lower(whichMAX)) %>% pull(id)
  b = tempocol %>% arrange(str_to_lower(x)) %>% pull(c)
  OMT <- givemereorderedMX(assay(dsv), a , b)
  
  mysplitcol = data.frame(x=sapply(colnames(OMT),
                                   function(x) unlist(str_split(x, "\\."))[2]))
  rownames(mysplitcol) = colnames(OMT)
  # pick  viridis scale for continuous values (vst)
  scalecontinuous = viridis(10, alpha = 1, begin = 1, end = 0, direction = 1)
  minscalecol = scalecontinuous[5]
  maxscalecol = scalecontinuous[1] 
  #col_fun = colorRamp2(c(0,20), c("lightgray", "firebrick4"))
  col_fun = colorRamp2(c(0,round(max(as.array(OMT))+1, 0)),
                       c(minscalecol, maxscalecol))
  agesannot = sapply(colnames(OMT), 
                     function(x) unlist(str_split(x, "\\."))[1])
  #listofcolorsage = sapply(agesvec, function(x) ifelse(x == "Old", "gold","darkgray"))
  hi = HeatmapAnnotation(age = agesannot,
                         col = list(age = c("Old"="gray10","Young"="gold")),
                         simple_anno_size = unit(3, "mm"))
  
  pdf(destinyfilepdf, width=6)
  labelsROW = Z[match(rownames(OMT), Z$id),]$whichMAX
  cchm = ComplexHeatmap::Heatmap(OMT,  col = col_fun,
                                 cluster_rows = F,
                                 cluster_columns = T,
                                 column_split = mysplitcol,
                                 row_split = labelsROW,
                                 row_names_gp = gpar(fontsize = 0 ),
                                 column_names_gp = gpar(fontsize = 0),
                                 column_names_rot = 45,
                                 show_row_names = F, 
                                 heatmap_legend_param = list(title="Transformed counts (vst)",
                                                             direction = "vertical"),
                                 bottom_annotation = hi
  )
  
   draw(cchm, 
         column_title = paste("Expression of genes with Tau >",
                              CUTOFFTAU, "across libraries (", d,")" ),
  ) 
  dev.off()
  
}

if (dopreplots){
  preliminaryheatmapvst('D0', paste0(odir, gdir, "1prep_Intx_heatmap_",'D0',".pdf") )
  preliminaryheatmapvst('D2', paste0(odir, gdir, "1prep_Intx_heatmap_",'D2',".pdf") )
  preliminaryheatmapvst('D4', paste0(odir, gdir, "1prep_Intx_heatmap_",'D4',".pdf") )
  preliminaryheatmapvst('D7', paste0(odir, gdir, "1prep_Intx_heatmap_",'D7',".pdf") )
}  
# ====================================== END preliminary heatmaps =====================================

# ==================== Test diff expr on Tau filtered  (if doDEtest is TRUE) ==================
print("explore multiple occurrences of ensembleids by day ")
for (d in daysv){
  print(d)
  foo <- table(dedup[[d]]$id)>=2
  print(paste(foo[foo==TRUE], "==> duplicated ensemblids"))
}
if (doDEtest){
  resEnsDE <- list()
  for (d in daysv){
    Z <- left_join(dedup[[d]], genes_df, by=c("id"="Geneid")) # ensemblids
    Z$symbol <- Z$symbol.x    # one symbol can correspond to two distinct ids 
    Z$symbol.y <- Z$symbol.x <- NULL
    smx <- fmat[rownames(fmat) %in% Z$id, str_detect(colnames(fmat),d)]
    smm <- metadata %>% filter(time==d)
    if (length(unique(rownames(smx))) == length(rownames(smx))){
      ( "ok" )
      keep.s <- apply(smx, 1, function(w) ifelse(sum(w >= 5)>= 3, T, F ))
      smx <- smx[keep.s,]
      ds.s <- DESeqDataSetFromMatrix(countData=smx, colData=smm, design = ~age)
      colData(ds.s)$age <- factor(colData(ds.s)$age,
                                  levels=c("Young","Old"))
      res <- DESeq(ds.s, full = ~age)
      restab <- as_tibble(results(res))
      restab$id = rownames(results(res))
      restab$symbol = Z[match(restab$id,Z$id),]$symbol
      resEnsDE[[d]] <- restab %>% filter(! is.na(padj)) 
      }else{
        print("stop, rownames must be set unique !!")
        stop()
      }
  }
  saveRDS(resEnsDE, paste0(odir, "rds/", "Intx_shot_onTauExtract.rds"))
}

# =========================== end Test diff Tau Filtered =======================

# ======================== GSEA on diff (gene symbols) =========================
resEnsDE <- readRDS(paste0(odir, "rds/", "Intx_shot_onTauExtract.rds"))

# pick genes of interest from resEnsDE full results list : 
infoinput <- data.frame("day" = c(), "inputsize"=as.integer(c()), 
                        "maxpadj" = c(), "minabslfc" = c())
intere_ <- list()
for (d in daysv){
  resUniqMx <- list()
  resUniqMx[[d]] <- resEnsDE[[d]] %>% group_by(symbol) %>% slice_max(baseMean, n=1)
  Z <- left_join(dedup[[d]], genes_df, by="symbol")
  resUniqMx[[d]]$celltype <- Z[match(resUniqMx[[d]]$symbol, Z$symbol),]$whichMAX
  intere_uplfc <- resUniqMx[[d]] %>% filter(log2FoldChange > 0) %>%
    group_by(celltype) %>% slice_max(log2FoldChange, n=M, with_ties=F)  %>%
    arrange(desc(log2FoldChange)) %>% filter(abs(log2FoldChange) >= 0.2)
  intere_downlfc <- resUniqMx[[d]] %>% filter(log2FoldChange < 0) %>%
    group_by(celltype) %>% slice_min(log2FoldChange, n=M, with_ties=F) %>%
    arrange(desc(log2FoldChange)) %>% filter(abs(log2FoldChange) >= 0.2)
  intere_[[d]] <- rbind(intere_uplfc, intere_downlfc)
  infoinput <- rbind(infoinput, c(d,
                                  dim(intere_[[d]])[1],
                                  max(intere_[[d]]$padj),
                                  min(abs(intere_[[d]]$log2FoldChange)) ) )
}
colnames(infoinput) <- c("day_celltype" , "inputsize", 
                         "maxpadj" , "minabslfc" )
infoinput$inputsize <- as.integer(infoinput$inputsize)
write.table(infoinput, paste0(odir,gdir, "csv/infoinput_Intx.csv"), sep='\t', 
            col.names = T, row.names = F)
dobarplotpregsea <- F
if (dobarplotpregsea){
    pdf(paste0(odir, gdir, "1prep_barplot_pregsea.pdf"), height = 8, width=12)
    for (d in daysv){
      gseagenes <- intere_[[d]] %>% pull(log2FoldChange)
      names(gseagenes) <- intere_[[d]]$symbol
      barplot(gseagenes, main = d)
    }
    dev.off()
}
# prepare gmt files for gsea 

print("Running GSEA on selected gene set")
Reac_gmt <- msigdbr(species = "Mus musculus", category = 'C2', subcategory=c('CP:REACTOME'))
msigdbr_list = split(x = Reac_gmt$gene_symbol, f = Reac_gmt$gs_name)
if (doGSEAintx){
  set.seed(42)
  fgseaR_ <- list()
  for (d in daysv){
    print(d)
    gseagenes <- intere_[[d]] %>% pull(log2FoldChange)
    names(gseagenes) <- intere_[[d]]$symbol
    fgseaRes <- fgsea::fgsea(pathways = msigdbr_list,
                             stats = gseagenes,
                             minSize = 5,
                             eps = 0,
                             maxSize = Inf
    ) # no nperm because Multilevel 
    fgseaR_[[d]] <- fgseaRes
  } 
}

rm(fgseaRes)
saveRDS(fgseaR_, paste0(odir, gdir, "rds/Intx_fgsea_full.rds"))




###  CAUTION:
##############################################################################
# this piece of code demonstrates how TPM and Counts can be discordant,
# let's take Vegfa gene at  M2 D2 :
prefil_cou <- "data/prefiltered_counts.rds"
fmat <- readRDS(prefil_cou)
foo <- sapply(colnames(fmat), function(x) str_detect(x, "M2.D2"))
vegfainfo <- list()
vegfainfo[["counts"]] <- fmat["ENSMUSG00000023951", foo]
vegfainfo[["age"]] <- sapply(names(foo[foo==T]), function(x) unlist(str_split(x, "\\."))[1])
ggplot(as_tibble(vegfainfo), aes(age,counts) ) +
  geom_point(aes(age,counts)) + 
  stat_summary( geom = "point", fun = "mean", col="white", size=4, shape=21,
                fill="salmon") +
  ggtitle("Counts, gene ENSMUSG00000023951 ( Vegfa ) in M2 D2")

prefil_tpm <- "data/prefiltered_TPM.rds"
fTPM <- readRDS(prefil_tpm)
vgf <- list()
vgf[["tpm"]] <- fTPM["ENSMUSG00000023951", foo]
vgf[["age"]] <- sapply(names(foo[foo==T]), function(x) unlist(str_split(x, "\\."))[1])
ggplot(as_tibble(vgf), aes(age,tpm)) + 
  geom_point(aes(age,tpm)) +
  stat_summary( geom = "point", fun = "mean", col="white", size=4, shape=21,
                fill="salmon") +
  ggtitle("TPM, gene ENSMUSG00000023951 ( Vegfa ) in M2 D2")
# THIS IS A UNRESOLVED QUESTION !!!! did not find biblio about it.
##############################################################################
#??end caution

