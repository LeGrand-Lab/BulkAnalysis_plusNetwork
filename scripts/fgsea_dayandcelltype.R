#  Runs GSEA by day and celltype
# from 'static' DEGs
# Saves rds objects into exam_INTER_conditions/static/GSEA/
# --
# johaGL
library(fgsea)
setwd("~/BulkAnalysis_plusNetwork/")
# take top 200 ranked genes by pvalue and abs(l2fc))
# and perform gsea like analysis with gprofiler2

odir = "exam_INTER_conditions/static/"
degfile = "shot_dataframe_softfilter.csv"   # this is soft filtered version

DEdf <- read.table(paste0(odir, degfile), header=T, sep='\t')
# DEfi <- read.table(paste0(odir, degfile), header=T, sep=',')
summary(DEdf)
# gene symbol missing, add it
genes_df <- read.table("data/genesinfo.csv",sep="\t",header=T)

DEdf$symbol <- genes_df[match(DEdf$id, genes_df$Geneid),]$symbol
head(DEdf)

# ======= split by day and keep in a list of dataframes (each df ~100 genes): =====
DE_l = list()  #  
for (k in c('D0','D2', 'D4', 'D7')){
  tmp <- DEdf %>% filter(day == k)  %>% 
    mutate(absLFC = abs(log2FoldChange)) %>% mutate(sens= ifelse(log2FoldChange < 0,"down", "up"))
  # keep top genes by abslfc and pval
  tmp_up <- tmp %>% filter(sens == "up") %>% group_by(type) %>%
    arrange(padj, desc(absLFC), .by_group = TRUE) %>% slice_min(padj, n=50)
  tmp_down <- tmp %>% filter(sens == "down") %>% group_by(type)   %>%
    arrange(padj, desc(absLFC), .by_group = TRUE) %>% slice_min(padj, n=50)
  DE_l[[k]] <-  rbind(tmp_up,tmp_down)
}

# ======================== perform Gsea on top ranked genes ====================
# the genes from previous step
thegmt <- read.table(paste0("stock_gmtfiles/","Hallmark_React.gmt"), sep='\t',
                     header=T)
msigdbr_list = split(x = thegmt$gene_symbol, f = thegmt$gs_name)

set.seed(42)

pathsFull_l <- list("D0"=list(),"D2"=list(),"D4"=list(),"D7"=list())
pathsFiltered  <- list("D0"=list(),"D2"=list(),"D4"=list(),"D7"=list())
for (k in c('D0', 'D2', 'D4', 'D7')){
  cts <- unique(DE_l[[k]]$type) 
  dd <- list()
  dd_f <- list()
  for (CT in cts){
    here.df <- DE_l[[k]] %>% filter(type == CT)
    print(c(k, CT))
    gseagenes = here.df %>% arrange(desc(absLFC)) %>% pull(log2FoldChange)
    names(gseagenes) <- here.df$symbol
    print(gseagenes)
    barplot(sort(gseagenes))
    pathsFull_l[[k]] = list()
    pathsFiltered[[k]] = list()
    fgseaRes <- fgsea::fgsea(pathways = msigdbr_list, 
                             stats = gseagenes,
                             minSize=3,
                             maxSize=Inf, nperm = 100000) 
    topPathwaysUp <- fgseaRes[ES>0][head(order(padj),n=15), ]
    topPathwaysDown <- fgseaRes[ES<0][head(order(padj), n=15), ]
    combipath <- rbind(topPathwaysUp, topPathwaysDown)
    combipath <- combipath %>% 
      mutate(sens = ifelse(ES > 0, "UP", "DOWN" ))
    dd[[CT]] <- fgseaRes
    dd_f[[CT]] <- combipath
  }
  pathsFull_l[[k]] <- dd
  pathsFiltered[[k]] <- dd_f
}

saveRDS(pathsFull_l, file=paste0(odir,"GSEA/fgseaByDay_full.rds" ))
saveRDS(pathsFiltered, file=paste0(odir, "GSEA/fgseaByDay_filtered.rds" ))


# ## initial test fgsea:
# set.seed(42)
# pathsFull_l = list()
# pathsFiltered = list()
# k = 'D0'
# cts <- unique(DE_l[[k]]$type) 
# CT <- 'ECs'
# here.df <- DE_l[[k]] %>% filter(type == CT)
# gseagenes = here.df %>% arrange(desc(absLFC)) %>% pull(log2FoldChange)
# names(gseagenes) <- here.df$symbol
# barplot(sort(gseagenes))
# pathsFull_l[[k]] = list()
# pathsFiltered[[k]] = list()
# fgseaRes <- fgsea::fgsea(pathways = msigdbr_list, 
#                          stats = gseagenes,
#                          minSize=3,
#                          maxSize=Inf, nperm = 100000) 
# 
# topPathwaysUp <- fgseaRes[ES>0][head(order(padj),n=15), ]
# topPathwaysDown <- fgseaRes[ES<0][head(order(padj), n=15), ]
# combipath <- rbind(topPathwaysUp, topPathwaysDown)
# combipath <- combipath %>% 
#   mutate(sens = ifelse(ES > 0, "UP", "DOWN" ))
# pathsFull_l[[k]][[CT]] <- fgseaRes
# pathsFiltered[[k]][[CT]] <- combipath
