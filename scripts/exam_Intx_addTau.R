# test to cross specificity and DE results
# by day by cell type
# JohaGL
setwd("~/BulkAnalysis_plusNetwork/")
library(tidyverse)
library(ggmosaic)

odir = "exam_INTER_conditions/static/"
DEf = readRDS(paste0(odir,"shot_rds_full.rds"))
dim(DEf)
head(DEf,2)

needconsensus <- F
doconcat <- F
daysv = c("D0","D2","D4","D7")
CUTOFFTAU = 0.5

# ============= creating consensus TAu lists (better use Extended version) =====
if (needconsensus){ 
  print("an extended version of consensus, including all Tau values")
  consensus_tau <- list()
  daysv = c("D0","D2","D4","D7")
  for (day in daysv){
    ty <- read.table(paste0("Tau/TauSpecificity_Young",day, ".txt"), sep="\t", header=T)
    to <- read.table(paste0("Tau/TauSpecificity_Old",day, ".txt"), sep="\t", header=T)
    # remember = Tau is calculated on meanTPM (mean over replicates' TPMs)
    #  therefore, its value is influenced by changes in expression due to age
    # for example, a hypothetical gene1 : [FAPS, sCs, ECs, M1], their meanTPMs
    # being :   Old =  [100,100,100,100]  ; Young = [100,2000,100,100]
    #  the result will say that gene1 is 'sCs' specific, but only in Young.
    #  As a matter of simplicity, by day and by cell type:
    #             if gene has max TPM in 1 single celltype in at least one age,
    # it will be included in consensus list:
    brut <- rbind(ty,to)
    brut <- brut  %>% filter(nbMAX == 1 ) %>%
      group_by(symbol, whichMAX) %>%    # whichMAX is the cell type
      slice(which.max(Tau)) 
    
    consensus_tau[[day]] <- brut %>% select(symbol, Tau, class, whichMAX)
  }
  for (i in (names(consensus_tau))){
    print("")
    print(i)
    print(table(consensus_tau[[i]]$whichMAX))
    print(table(consensus_tau[[i]]$class))
    
  }
  saveRDS(consensus_tau, paste0(odir,"conseTau4DEGs/", "Ext_conseTau.rds"))
}else{
  print(paste0("consensuslist already exists in :",
               odir,"conseTau4DEGs/", "Ext_conseTau.rds"))}
# ===========================end creating consensu =============================


# ============================ give uniqueness accross =============================
# UNIQUE list:
consensus_tau <- readRDS(paste0(odir,"conseTau4DEGs/", "Ext_conseTau.rds"))
# extract only Tau > 0.3 and deduplicate with max tau
dedup <- list()
for (d in daysv){
  tf <- consensus_tau[[d]] %>% filter(Tau >= 0.3) %>% 
      group_by(symbol) %>% slice_max(Tau)
  dedup[[d]] <- tf
}

# # bring counts matrix, and by day, do pca only on 
# respca <- list()
#   dedup[[d]]$symbol 
# }

fmat <- readRDS("data/prefiltered_counts.rds")
metadata <- readRDS("data/metadata.rds")
genes_df <- read.table("data/genesinfo.csv", sep="\t", header=T)

d = "D2"
h <- left_join(dedup[[d]], genes_df, by="symbol")
names(h) <- dedup[[d]]$whichMAX

smx <- fmat[rownames(fmat) %in% h$Geneid, str_detect(colnames(fmat),d)]
smm <- metadata %>% filter(time==d)
library(DESeq2)
ds <- DESeqDataSetFromMatrix(countData=smx, colData=smm, design = ~age)
dsv <- DESeq2::varianceStabilizingTransformation(ds)
library(pheatmap)
pheatmap(assay(dsv), cluster_c)
library(ComplexHeatmap)

mysplitcol = data.frame(x=sapply(colnames(assay(dsv)), 
                                 function(x) unlist(str_split(x, "\\."))[2]))
rownames(mysplitcol) = colnames(assay(dsv))
library(circlize)
col_fun = colorRamp2(c(0,20), c("steelblue", "firebrick4"))
pdf(paste0(odir, "VSTtaufiltered.pdf"), width=15 )

agesvec = sapply(colnames(assay(dsv)), 
                 function(x) unlist(str_split(x, "\\."))[1])
listofcolorsage = sapply(agesvec, function(x) ifelse(x == "Old", "gold","darkgray"))

ho = HeatmapAnnotation(foo = anno_block(gp = gpar(fill = listofcolorsage )))
cchm = ComplexHeatmap::Heatmap(assay(dsv),  col = col_fun,
                              cluster_rows = TRUE,
                              cluster_columns = TRUE,
                              bottom_annotation = ho,
                              #column_split = mysplitcol,
                              row_names_gp = gpar(fontsize = 0 ),
                              column_names_gp = gpar(fontsize = 9),
                              column_names_rot = 45,
                              heatmap_legend_param = list(title="Expression (vst)",
                                                          direction = "vertical")
                             
                          )

draw (cchm, column_title = "Libraries and genes with Tau > 0.3" )

dev.off()

## todo : verify colors old, colors young (gold vs black)
## and then verify annotation and legend
# ================================================================================
      

# ================ concatenate Tau and *full* DEG information =========================
if (doconcat){
  #consensus_tau <- readRDS( paste0(odir,"conseTau4DEGs/", "conseTau.rds"))
  consensus_tau <- readRDS( paste0(odir,"conseTau4DEGs/", "Ext_conseTau.rds") )
  # cross information with DEG results FULL version
  degfull <- readRDS(paste0(odir, "shot_rds_full.rds"))
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
  saveRDS(giant, paste0(odir, "TauPlusDEinfo_full.rds"))
}

# ===============================  load element for plots:   :=======================
giant <- readRDS(paste0(odir, "TauPlusDEinfo_full.rds"))
for (d in daysv){
  f <- giant %>% filter(day==d & Tau_discrete=="SPECIFIC") 
  n_occ <- data.frame(table(f$symbol))
  print("checking each day , printing duplicates")
  print(unique(f$type))
  if(dim(n_occ[n_occ$Freq > 1,])[1] == 0){ print("ok for plots (no duplics)")
    }else{print("attention, duplicated symbols across celltypes this day")}
}
# =================================  end load =======================

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



