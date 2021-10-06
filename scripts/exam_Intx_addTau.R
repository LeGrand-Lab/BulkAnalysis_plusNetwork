# Test Differential enrichment of pathways
# BY DAY, on unique ensemblids based on Tau filtering
# outputs into 'exam_inter_conditions/static' :
#           rds/Intx_conseTau_ensemblid.rds
#           rds/Intx_shot_onTauExtract.rds
#           GSEA/rds/Intx_fgsea_full.rds
#           
# JohaGL
setwd("~/BulkAnalysis_plusNetwork/")
library(tidyverse)
library(ggmosaic)
library(DESeq2)

odir = "exam_INTER_conditions/static/"
taudir = "Tau/"
consensusfile <- "conseTau_ensemblid.rds"

#taudefile <- "TauPlusDEinfo_full.rds"
taudefile <- "TauPlusDEinfo_fullV4.rds"
doDEtest <- F
daysv = c("D0","D2","D4","D7")
CUTOFFTAU = 0.3
fmat <- readRDS("data/prefiltered_counts.rds")
metadata <- readRDS("data/metadata.rds")
genes_df <- read.table("data/genesinfo.csv", sep="\t", header=T)



# ======== Import tauconsensus and give uniqueness: function importdedup() =======
importdedup <- function(consensusfile, CUTOFFTAU){
    consensus_tau <- readRDS(paste0(taudir, consensusfile))
    # extract only Tau > CUTOFFTAU and deduplicate with max tau
    dedup <- list()
    for (d in daysv){
      tf <- consensus_tau[[d]] %>% filter(Tau >= CUTOFFTAU) %>% 
        group_by(id) %>% slice_max(Tau)  ### ??? is this ok ????
      dedup[[d]] <- tf
      } # end for
    return(dedup)
  }
# =============================================================================

# ===================== heatmap day by day ====================================
library(ComplexHeatmap)

library(circlize)
library(viridis)
# !!  use here function importdedup
dedup <- importdedup(consensusfile, CUTOFFTAU)


d = "D7" ############ !!!!!!!!!!!!!!!!!!!!!!! by day !!!!!!!!!!!!!!!!!
zo <- left_join(dedup[[d]], genes_df, by=c("id"="Geneid")) # ensemblids
zo$symbol <- zo$symbol.x    # one symbol can correspond to two distinct ids 
zo$symbol.y <- zo$symbol.x <- NULL
smx <- fmat[rownames(fmat) %in% zo$id, str_detect(colnames(fmat),d)]
smm <- metadata %>% filter(time==d)
#  rownames for this day must be unique 
if (length(unique(rownames(smx))) == length(rownames(smx))){
  print("ok, rownames are unique (ensemblid)")
}else{
  print("stop, rownames must be set unique")
}

predsv <- DESeqDataSetFromMatrix(countData=smx, colData=smm, design = ~age)
dsv <- DESeq2::varianceStabilizingTransformation(predsv)
dsv <- DESeq2::rlog(predsv)
mysplitcol = data.frame(x=sapply(colnames(assay(dsv)),
                                 function(x) unlist(str_split(x, "\\."))[2]))
rownames(mysplitcol) = colnames(assay(dsv))
# reorder matrix, begin with a empty, then give rownames :

orderrows = zo[match(rownames(assay(dsv)), zo$id),]

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
FOFO <- givemereorderedMX(assay(dsv), a , b)

mysplitcol = data.frame(x=sapply(colnames(FOFO),
                                 function(x) unlist(str_split(x, "\\."))[2]))
rownames(mysplitcol) = colnames(FOFO)


# pick  viridis scale for continuous values (vst)
scalecontinuous = viridis(10, alpha = 1, begin = 1, end = 0, direction = 1)
minscalecol = scalecontinuous[5]
maxscalecol = scalecontinuous[1] 
#col_fun = colorRamp2(c(0,20), c("lightgray", "firebrick4"))
col_fun = colorRamp2(c(0,round(max(as.array(FOFO))+1, 0)),
                     c(minscalecol, maxscalecol))
agesannot = sapply(colnames(FOFO), 
                 function(x) unlist(str_split(x, "\\."))[1])
#listofcolorsage = sapply(agesvec, function(x) ifelse(x == "Old", "gold","darkgray"))
hi = HeatmapAnnotation(age = agesannot,
                       col = list(age = c("Old"="gray10","Young"="gold")),
                       simple_anno_size = unit(3, "mm"))

pdf(paste0(odir, "exam_Intx_heatmap_",d,".pdf"), width=6)
KKK = zo[match(rownames(FOFO), zo$id),]$whichMAX
cchm = ComplexHeatmap::Heatmap(FOFO,  col = col_fun,
                              cluster_rows = F,
                              cluster_columns = T,
                              column_split = mysplitcol,
                             # split = KKK,
                              row_split = KKK,
                              row_names_gp = gpar(fontsize = 0 ),
                              column_names_gp = gpar(fontsize = 0),
                              column_names_rot = 45,
                              show_row_names = F, 
                              heatmap_legend_param = list(title="Transformed counts (vst)",
                                                          direction = "vertical"),
                              bottom_annotation = hi
                          )

draw ( cchm, 
      column_title = paste("Expression of genes with Tau >",
                           CUTOFFTAU, "across libraries (", d,")" ),
      )
dev.off()

## and then verify annotation and legend
# ============================ END heatmap =====================================

# ==================== Test diff expr on Tau filtered  =========================
dedup <- importdedup(consensusfile, CUTOFFTAU)
for (d in daysv){
  print(d)
  foo <- table(dedup[[d]]$id)>=2
  print(paste(foo[foo==TRUE], "==> duplicated ensemblids"))
}
if (doDEtest){
  resEnsDE <- list()
  for (d in daysv){
    zo <- left_join(dedup[[d]], genes_df, by=c("id"="Geneid")) # ensemblids
    zo$symbol <- zo$symbol.x    # one symbol can correspond to two distinct ids 
    zo$symbol.y <- zo$symbol.x <- NULL
    smx <- fmat[rownames(fmat) %in% zo$id, str_detect(colnames(fmat),d)]
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
      restab$symbol = zo[match(restab$id,zo$id),]$symbol
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

library(fgsea)
# pick genes of interest from resEnsDE full results list : 
infoinput <- data.frame("day" = c(), "inputsize"=as.integer(c()), 
                        "maxpadj" = c(), "minabslfc" = c())
intere_ <- list()
for (d in daysv){
  # !! attention: there can be symbols repetitions (as DE ran on ensemblids)
  # make unique (slice_max on baseMean): 
  resUniqMx <- list()
  resUniqMx[[d]] <- resEnsDE[[d]] %>% group_by(symbol) %>% slice_max(baseMean, n=1)
  zo <- left_join(dedup[[d]], genes_df, by="symbol")
  resUniqMx[[d]]$celltype <- zo[match(resUniqMx[[d]]$symbol, zo$symbol),]$whichMAX
  intere_uplfc <- resUniqMx[[d]] %>% filter(log2FoldChange > 0) %>%
    group_by(celltype) %>% slice_max(log2FoldChange, n=1000, with_ties=F)  %>%
    arrange(desc(log2FoldChange)) %>% filter(abs(log2FoldChange) >= 0.2)
  intere_downlfc <- resUniqMx[[d]] %>% filter(log2FoldChange < 0) %>%
    group_by(celltype) %>% slice_min(log2FoldChange, n=1000, with_ties=F) %>%
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
write.table(infoinput, paste0(odir,"GSEA/infoinput_Intx.csv"), sep='\t', 
            col.names = T, row.names = F)
dobarplotpregsea <- F
if (dobarplotpregsea){
    pdf(paste0(odir, "GSEA/intx_barplot_pregsea.pdf"), height = 8, width=12)
    for (d in daysv){
      gseagenes <- intere_[[d]] %>% pull(log2FoldChange)
      names(gseagenes) <- intere_[[d]]$symbol
      barplot(gseagenes, main = d)
    }
    dev.off()
}
# prepare gmt files for gsea 
library(msigdbr)
Reac_gmt <- msigdbr(species = "Mus musculus", category = 'C2', subcategory=c('CP:REACTOME'))
msigdbr_list = split(x = Reac_gmt$gene_symbol, f = Reac_gmt$gs_name)
doGSEAintx <- F
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
                             maxSize = Inf
    ) # no nperm because Multilevel 
    fgseaR_[[d]] <- fgseaRes
  } 
}

rm(fgseaRes)
saveRDS(fgseaR_, paste0(odir, "GSEA/rds/Intx_fgsea_fullNEW.rds"))
# plots 
for (d in daysv){
  # pending to fillllll
}

library(gridExtra)
library(cowplot)

# DO PLOTS, ORGANIZE BY PADJ (from smallest)
for (d in daysv){
  pdf(paste(odir, "meh", d,".pdf" ), width=15, height=15)
  fgseaRes <- fgseaR_[[d]]
  summary(fgseaRes$padj)
  NumP = 15
  topPathwaysUp <- fgseaRes[ES>0][head(order(padj),n=NumP), ]
  topPathwaysDown <- fgseaRes[ES<0][head(order(padj), n=NumP), ]
  topPathBoth <- c(topPathwaysUp$pathway, rev(topPathwaysDown$pathway))
  ouif = fgsea::plotGseaTable(msigdbr_list[topPathBoth], gseagenes, fgseaRes, 
                              gseaParam = 0.5 , render=F) 

  plot_grid(ggdraw() + draw_label(paste(d , ": Top enriched Pathways (GSEA), Old vs Young")),
            plot_grid(NULL, ouif , rel_widths =c(4,7)),
            NA,
            nrow = 3, rel_heights = c(1,9,9))
  plotsenrichu_ <- list()
  for (i in 1:NumP){
    pup = topPathwaysUp[i,]
    plotsenrichu_[[i]] <- fgsea::plotEnrichment(msigdbr_list[[pup$pathway]], gseagenes) + 
      labs(title= d, subtitle = pup$pathway, 
           caption=paste( "(NES:", pup$NES, ", padj :", pup$padj,")"))
  }
  plot_grid(plotlist = plotsenrichu_, nrow = 5)
  plotsenrichdw_ <- list()
  for (i in 1:NumP){
    pdw = topPathwaysDown[i,]
    plotsenrichdw_[[i]] <- fgsea::plotEnrichment(msigdbr_list[[pdw$pathway]], gseagenes) + 
      labs(title= d, subtitle = pdw$pathway,
           caption=paste( "(NES:", pdw$NES, ", padj :", pdw$padj,")") )
  }
  plot_grid(plotlist = plotsenrichdw_, nrow = 5)
  # dos <- fgsea::plotEnrichment(msigdbr_list[[topPathBoth[1]]], gseagenes) + 
  #   labs(title= d, subtitle=topPathwaysUp[1,]$pathway,
  #        caption=paste(d,"Old vs Young" ))
  # tres <- fgsea::plotEnrichment(msigdbr_list[[topPathBoth[NumP*2]]], gseagenes) + 
  #   labs(title= d, subtitle=topPathwaysDown[1,]$pathway,
  #        caption=paste(d,"Old vs Young"))
  # plot_grid(dos, tres, NA, ncol=3, rel_widths = c(4,4,2))
  dev.off()
}


# ==============================================================================
      

# ================ concatenate Tau and *full* DEG information ==================
# CONSIDER TO TRANSFORM THIS THIS PART INTO MARKDOWN
# This uses DEG that resulted from script exam_Inter_cond_sta.R
if (doconcat){
  consensus_tau <- readRDS( paste0(odir,"rds/", consensusfile) )
  # cross information with DEG results FULL version
  degfull <- readRDS(paste0(odir, "rds/shot_rds_full.rds"))
  daysv = c("D0","D2","D4","D7")
  mdegtau <- list()
  for (d in daysv){
  D0deg <- degfull %>% filter(day == d) %>% select(symbol,log2FoldChange,padj,type)
  
  df_tmp <- right_join(consensus_tau[[d]] %>% mutate(type = whichMAX),
                      D0deg %>% filter(padj <= 0.5), # this also excludes NA padj
                      by=c("symbol","type"))
  
  df_tmp <- df_tmp %>% filter(abs(log2FoldChange) > 0.1)
  df_tmp <- df_tmp %>% mutate(tauAvailable=ifelse(is.na(Tau), "no", "yes"))
  df_tmp <- df_tmp %>% mutate(Tau_discrete = case_when(
   # is.na(Tau) ~ "NA", # was not available
        Tau >= CUTOFFTAU ~ "SPECIFIC",
        TRUE ~ "HOUSEKEEPING"
      ), LFC_discrete = case_when(
        log2FoldChange >= 1 ~ "A. UP",
        log2FoldChange < 1 & log2FoldChange >= .5  ~ "B. 0.5 to 1",
        log2FoldChange <= -.5 & log2FoldChange > -1 ~ "D. -1 to -0.5",
        log2FoldChange <= -1 ~  "E. DOWN",
        TRUE ~ "C. -0.5 < lfc < 0.5"     )  ) # end mutate
  
  df_tmp <- df_tmp %>% mutate(padj_discrete =  case_when(
    is.na(padj) ~ "NApadj",
    padj <= 0.005 ~ "a. SIGNIFICANT (<= 0.005)",
    padj <= 0.05 & padj > 0.005 ~ "b. SIGNIFICANT (0.05 to 0.005)",
    padj > 0.05 & padj <= 0.5 ~ "c. 0.05 to 0.5",
    padj > 0.5 ~ "d. more than 0.5") )
  df_tmp$day <- d
  mdegtau[[d]] <- df_tmp
  }
  giant = bind_rows(mdegtau)
  # check EACH DAY: there are still some genes classified as "specific" that appear 
  # *duplicated*  among two cell types
  for (d in daysv){
    f <- giant %>% filter(day==d & Tau_discrete=="SPECIFIC") 
    n_occ <- data.frame(table(f$symbol))
    print("checking each day , printing duplicates")
    print(unique(f$type))
    print(n_occ[n_occ$Freq > 1,])
  }
  # make unique (a symbol cannot be repeated across celltypes in one day )
  # criterion : smallest padj
  mdegtau <- list() # clear to reuse this var
  for (d in daysv){
    f <- giant %>% filter(day==d) %>% group_by(symbol) %>%
      slice_min(padj)
    mdegtau[[d]] <- f
  }
  giant <- bind_rows(mdegtau)
  saveRDS(giant, paste0(odir, "rds/", taudefile))
}

# ===============================  load element for plots:   :=======================
giant <- readRDS(paste0(odir, "rds/", taudefile))
for (d in daysv){
  f <- giant %>% filter(day==d & Tau_discrete=="SPECIFIC") 
  n_occ <- data.frame(table(f$symbol))
  print("checking each day , printing duplicates")
  print(unique(f$type))
  if(dim(n_occ[n_occ$Freq > 1,])[1] == 0){ print("ok for plots (no duplics)")
    }else{print("attention, duplicated symbols across celltypes this day")}
}
# =====================================  end load ==============================

# ==================================== MAplot ===============================
giant <- readRDS(paste0(odir, "rds/", taudefile))
#pdf(paste0(odir, "TauPlusDEinfo_plot.pdf"), width=10, height=6)
ggplot(giant, aes(x=Tau, y=log2FoldChange, color= Tau >= .5) ) + 
  geom_point(size=.5) + geom_hline(yintercept = c(-1.2,1.2),
                                   linetype = "dashed") +
  scale_color_brewer(palette = "Dark2") + facet_grid(padj<=0.05 ~ day) + 
  labs(title="DE information and Tau index", 
       caption="y axis plot labels: TRUE correspond to padj <=0.05") 
#dev.off()
# ====================================  end MAplot ===========================


# =================================  use Mosaic Plot :=======================
# reorder factors : 
giant$Tau_discrete <- ordered(giant$Tau_discrete, levels = c("SPECIFIC", "HOUSEKEEPING"))
revlevelpj <- rev(levels(as.factor(giant$padj_discrete) ))
giant$padj_discrete <- ordered(giant$padj_discrete, levels = revlevelpj)
#revlevelLFC <- rev(levels(as.factor(giant$LFC_discrete) ))
#giant$LFC_discrete <- ordered(giant$LFC_discrete, levels = revlevelLFC)
modx_labs = sapply(levels(as.factor(giant$padj_discrete)), function(x)
  unlist(str_split(x,"\\."))[1])

pdf(paste0(odir,"TauPlusDEinfo_mosaic_C.pdf"), width=12, height=10)
ggplot(giant)  +
geom_mosaic(aes(x=product(LFC_discrete), 
                conds = product(Tau_discrete), 
                fill = padj_discrete )) + 
  scale_fill_manual(values=c( "#404788FF", "#238A8DFF",
                                 "#55C667FF"), 
                    name="padj") +
  labs(x= "padj", title="Tau, Log2FC and padj values") + 
   scale_x_productlist("padj_categorical", labels=c("c","b","a", "","","")) +
   theme_mosaic() + theme(axis.text.x = element_text(angle = 90),
        plot.margin = unit(c(.1,.1,.1,.1), "cm"),
        plot.title= element_text(hjust=0.5)) + coord_flip() +  
  facet_grid(day~Tau_discrete)
dev.off()
  
# manual colors viridis --> scale_fill_manual(values=c("#440145FF", "#404788FF", "#238A8DFF",
#                              "#55C667FF", "#FDE725FF")),   
# =================================  another Plot :=======================
# proportion of housekeeping and prop specific along celltypes and days

worko <- giant %>% filter(abs(log2FoldChange) >= 0.5) # select only interesting changes

maxpj = round(max(worko$padj),1)
minlfc = round(min(abs(worko$log2FoldChange)),1)

ggplot(worko %>% group_by(type) %>% tally(),
       aes(x=type,y=n)) + geom_bar(position="dodge", stat="identity") + 
  geom_text(aes(label=n), position=position_dodge(width=0.9),vjust=-.25) +
  labs(title=paste("Number of top listed genes (padj < ",maxpj,", abs(lfc) ",minlfc,")"))

# among top listed genes (as seen in barplot)
# check proportion of "specific" Tau_discrete
ox <- worko %>% group_by(time_Day, type, Tau_discrete) %>% 
         summarise(n = n()) %>%
        group_by(time_Day,type) %>%
        mutate(Prop = n / sum(n))
ggplot(ox,
       aes(x=time_Day, y=Prop, fill=Tau_discrete)) +
  geom_col() + facet_grid(~type) +
labs(title = "Tissue Specific Proportion among Top listed genes, by cell type") 

# oo <- worko %>%
#   group_by(type, time_Day, class, padj_discrete, 
#                          LFC_discrete, Tau_discrete) %>%
#    summarise(sumTau = sum(Tau), 
#              sumAbsLFC = sum(abs(log2FoldChange)), 
#              sumpadj = sum(padj),
#              n=n())

# cts <- unique(giant$type)
# for (ct in cts){
#   for (d in  c("D0","D2","D4","D7")){
#     f <- giant %>% filter(type==ct, day==d)
#     print(ct)
#     print(length(unique(f$symbol))==length(f$symbol))
#   }
# }

# =========== test geom_mosaic
data(titanic)

ggplot(data= titanic) +
  geom_mosaic(aes(x = product(Class), 
                  conds = product(Age), fill = Survived)) +
  scale_fill_manual(values=c("#440145FF",  "#238A8DFF")) +
  theme_mosaic()


#  ========================================microtests complexheatmap=====

nr1 = 4; nr2 = 8; nr3 = 6; nr = nr1 + nr2 + nr3
nc1 = 6; nc2 = 8; nc3 = 10; nc = nc1 + nc2 + nc3
mat = cbind(rbind(matrix(rnorm(nr1*nc1, mean = 1,   sd = 0.5), nr = nr1),
                  matrix(rnorm(nr2*nc1, mean = 0,   sd = 0.5), nr = nr2),
                  matrix(rnorm(nr3*nc1, mean = 0,   sd = 0.5), nr = nr3)),
            rbind(matrix(rnorm(nr1*nc2, mean = 0,   sd = 0.5), nr = nr1),
                  matrix(rnorm(nr2*nc2, mean = 1,   sd = 0.5), nr = nr2),
                  matrix(rnorm(nr3*nc2, mean = 0,   sd = 0.5), nr = nr3)),
            rbind(matrix(rnorm(nr1*nc3, mean = 0.5, sd = 0.5), nr = nr1),
                  matrix(rnorm(nr2*nc3, mean = 0.5, sd = 0.5), nr = nr2),
                  matrix(rnorm(nr3*nc3, mean = 1,   sd = 0.5), nr = nr3))
)
mat = mat[sample(nr, nr), sample(nc, nc)] # random shuffle rows and columns
rownames(mat) = paste0("row", seq_len(nr))
colnames(mat) = paste0("column", seq_len(nc))

Heatmap(mat, name = "mat", col = col_fun, 
        column_title = "set a color vector for a continuous matrix")
#===================================================

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

##############################################################################
# end caution

# BIN:
# plotlfc <- 4
# MAplots[[d]] <- ggplot(data = restab %>% 
#                          mutate(log2FoldChange = case_when(
#                            log2FoldChange > plotlfc ~ plotlfc,
#                            log2FoldChange < -plotlfc ~ -plotlfc,
#                            TRUE ~ log2FoldChange
#                          ))) +   geom_point(
#   aes(x = baseMean, y = log2FoldChange, fill=padj),
#   size = .5,
#   #color = ifelse(restab$padj < 0.1, "red", "gray"),
#   alpha = 0.5    ) +
# scale_fill_continuous(type="viridis", ) +
# scale_x_log10()  

# --
# par(mfrow=c(2,2))
# for (d in daysv){
#   zo <- left_join(dedup[[d]], genes_df, by="symbol")
#   tmp <- DESeqDataSetFromMatrix(countData = fmat[rownames(fmat) %in% zo$Geneid, 
#                                                  str_detect(colnames(fmat),d)],
#                                 colData = metadata %>% filter(time==d),
#                                 design= ~age)
#   tmp$age <- relevel(tmp$age, ref="Young")
#   da <- estimateSizeFactors(tmp)
#   da <- estimateDispersions(da)
#   plotDispEsts(da, main=d)
#   rm(tmp)
# }
# par(mfrow=c(1,1))
## -- 
