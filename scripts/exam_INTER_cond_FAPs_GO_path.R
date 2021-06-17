# JohaGL 2021
# ---
library(dplyr)
library(tidyverse)
library(ggplot2)
library(cowplot)
library(RColorBrewer)
library(DESeq2)
library(ggsci) # publishing palettes
library(cowplot)
library(gridExtra)
library(reshape2)
library(ggthemes)
library(scales) # 'viridis_pal' et autres pal
library("BiocParallel")
register(MulticoreParam(12)) # TODO:  set n of cores depending of available
library(gprofiler2)
library(msigdbr)
library(fgsea)
setwd("~/BulkAnalysis_plusNetwork/")
# NOTE: color blindness for scales picked from:
#https://davidmathlogic.com/colorblind/#%23D81B60-%231E88E5-%23FFC107-%23004D40

ct = "FAPs"
gostneeded = F
go_path_doplots = T
gologfoldcutoff = 1.5
gopadjcutoff =  0.005
topgopath = 12

prefil_cou <- "data/prefiltered_counts.rds"
metadata.rds <- "data/metadata.rds"
genes_df <- read.table("data/genesinfo.csv",sep="\t",header=T)

resdir = "exam_INTER_conditions/dynamic/"
fmat <- readRDS(prefil_cou)
metadata <- readRDS(metadata.rds)
metadata <- metadata %>% mutate(typetimeage = paste0(type,".",time,".",age))
##  set which rows to keep :
keep <- apply(fmat, 1, function(row) ifelse(count(row >=5)>= 3, TRUE, FALSE) )
fmat <- fmat[keep,]

all_g_df <- read.table(paste0(resdir, ct, "_INTERagetime.csv"),sep='\t', header = T)

# =============== prepare data   =================
g_df <- all_g_df %>% filter(abs(log2FoldChange) >= gologfoldcutoff & 
                              padj <= gopadjcutoff ) 
initmeta = metadata %>% filter(type==ct) %>% select(-c(group,sample))
# ==================================================================

querygenes <- g_df$padj
names(querygenes) <- g_df$id
rankedgenes <- sort(querygenes)

# ==================  Gprofiler2 pathways enrichment: ========================
if (gostneeded){
  ctres = gprofiler2::gost(names(rankedgenes), organism='mmusculus',
                           ordered_query = TRUE, significant = TRUE,
                           user_threshold = 0.05, evcodes = TRUE,
                           correction_method = "fdr", 
                           sources = c('GO:BP', 'GO:MF', 'GO:CC' , 'KEGG', 
                                       'REAC', 'TF', 'MIRNA', 'CORUM', 'HP', 'HPA', 'WP') )
  View(ctres$result)
  saveRDS(ctres, paste0(resdir, ct,"_go_rawobject.rds")) #saved result
  
  # prepare for table
  # -----------------------------------------------------------------
  gp_mod = ctres$result %>% filter(precision > 0.04 & recall > 0.04) %>%
    select(query, source, term_id,
           term_name, p_value,query_size, 
           intersection_size,term_size,
           effective_domain_size,
           evidence_codes, intersection) 
  names(gp_mod) = c("Cluster", "Category", "ID", "Description", "FDR",
                    "query_size", "Count", "term_size", "effective_domain_size",
                    "evidence_codes", "genesEnriched" )
  
  print("all gprofiler2 significant results (databases), for this celltype: ")
  print(unique(gp_mod$Category))
  gp_mod = as_tibble(gp_mod) %>% 
    mutate(GeneRatio = Count/query_size,
           BgRatio = term_size/effective_domain_size,
           celltype=ct) %>% filter(FDR <= 0.05)
  #  exclude KEGG (kegg is too much redundant) , and Helicobacter Pylori  
  gp_mod = gp_mod %>%  filter(!Category %in% c("KEGG", "HP"))
  # 20 relevant pathways by category , to save into table : 
  gp_mod= gp_mod %>% group_by(Category)  %>%
    slice_min(order_by = FDR, n = 20) %>% 
    select(Category, ID, Description, FDR, genesEnriched,  
           query_size, Count, BgRatio, GeneRatio,  celltype ) 
  vec_genesymbols <- sapply(gp_mod$genesEnriched, function(x) {
    eids <- unname(unlist(str_split(x, ",")))
    symbs <- genes_df[match(eids, genes_df$Geneid),]$symbol
    return(unname(paste0(symbs, collapse=',')))
  })
  gp_mod$geneSymbols = vec_genesymbols
  write.table(gp_mod %>% relocate(geneSymbols, .after=FDR), 
              paste0(resdir, "go_gsea_csv/",ct, "_go_filtered.csv"), 
              sep='\t', col.names = T, row.names = F)
  # -------------------------------------------------------------------
  # end saved table
} # end if gostneeded
# NOTE:  last query on gprofiler2  : 17-06-2021 

# ======================== do plots go_paths ==================================
if (go_path_doplots){
  gp_mod = read.table(paste0(resdir, "go_gsea_csv/",ct, "_go_filtered.csv"),
                      header=T, sep='\t')
  GG <-  gp_mod %>% 
    mutate(enriched_terms = str_trunc(Description, 50, "right")) %>%
    group_by(Category) %>%
    slice_min(order_by = FDR, n = topgopath )
  # if very close terms, extract the one with max Count (using enriched_terms)
  GG <- GG %>% group_by(enriched_terms) %>% slice_max(Count, n=1, with_ties=F)
  gg_go <- ggplot(GG,   aes(enriched_terms, Count, fill=FDR)) +
    geom_bar(stat="identity", width=0.8) +
    scale_fill_continuous(type = "viridis", direction=-1, alpha=.8) +
    facet_grid(Category~.,scale="free", space="free") +
    coord_flip() +
    theme(axis.text.y = element_text(size=7)) +
    theme_bw() + theme( legend.title = element_text(size=7), 
                        legend.text  = element_text(size=7),
                        legend.key.size = unit(1, "lines")) +
    scale_y_continuous(breaks=breaks_pretty()) +
    labs(subtitle = paste("GO/Pathways enrichment"),
         caption = paste(ct, "Old vs Young, pathway enrichment significance (FDR)",
                         '\nabslfc >=', gologfoldcutoff, 
                         '  p<=',gopadjcutoff, '  n=',length(querygenes)) ) 
  # get genes gathered from all pathways
  toppathgenes = unique(unlist(str_split(GG$genesEnriched, ","))) 
  #pick genes to plot: genes having padj eq or inf to those found enriched:
  vizpadj = max(g_df[g_df$id %in% toppathgenes, ]$padj) # max padj enrich genes
  viznbtop = dim(g_df %>% filter(padj<=vizpadj))[1] + 2 
  VIZ <- g_df %>% slice_min(order_by = padj, n = viznbtop)
  VIZ$symbol <- genes_df[match(VIZ$id,genes_df$Geneid),]$symbol
  VIZ$wasenriched <- ifelse(VIZ$id %in% toppathgenes, 1, 0)
  
  # build labeller
  mylabeller <- paste0(VIZ$symbol," ", "(",
                       round(VIZ$log2FoldChange,1),"|",
                       round(VIZ$padj,2),")" )
  names(mylabeller) <- VIZ$symbol
  
  dso <- DESeqDataSetFromMatrix(fmat, metadata,
                                design= ~ age + time + age:time)
  subset <- dso[ unique(g_df$id), dso$type==ct ]
  subvsd <- varianceStabilizingTransformation(subset)
  zscores = t(scale(t(subvsd@assays@data@listData[[1]])))
  zscogene <- array(zscores[VIZ$id, ], 
                    dim= c(length(VIZ$id), dim(zscores)[2]))
  colnames(zscogene) <- colnames(zscores)
  print("inserting symbols instead  of ensembl for plotting")
  rownames(zscogene) <- VIZ$symbol
  if(all(rownames(initmeta)==colnames(zscogene))){
    print("rownames(initmeta)==colnames(zscores), ok for reshape2::melt")
    mt = cbind(initmeta,(t(zscogene)))
    mtmelt = reshape2::melt(mt, id.vars=c("newname", "time","type", "age",
                                          "typetimeage"),
                            variable.name="symbol", value.name="Zscore")  
    mtmelt$day = as.numeric(sapply(mtmelt$time, function(x)str_replace(x, "D","")))
    mtmelt$wasenriched = VIZ[match(mtmelt$symbol, VIZ$symbol),]$wasenriched
  }else{print("error, zcores colnames and rownames subsetted metadata NOT EQUAL")}
  
  gg_genenri <- ggplot(mtmelt %>% filter(wasenriched==1), 
                       aes(x=day, y=Zscore, group=age , color=age))  +
    geom_point(size=.4, aes(color=age)) + 
    geom_smooth(se=T,
                method="loess", 
                size=0.5,
                alpha=.1) + 
    scale_color_colorblind() + 
    scale_x_continuous(breaks=seq(min(mtmelt$day), max(mtmelt$day), 2)) +
    facet_wrap(~symbol, labeller = labeller(symbol=mylabeller),
               ncol = 5) + 
    theme_light() + theme(legend.position = c(1, 0),
                          legend.justification = c(1,0)) +
    labs(subtitle=paste(ct,"Old vs Young, GO/pathways retrieved"),
         caption=paste('total query (n genes):',dim(g_df)[1] ,
                       '\nGeneSymbol (logFC|padj)'))
  # genes not retrieved in enrichment for functional terms
  gg_genenon <- ggplot(mtmelt %>% filter(wasenriched==0), 
                       aes(x=day, y=Zscore, group=age , color=age))  +
    geom_point(size=.4, aes(color=age)) + 
    geom_smooth(se=T,
                method="loess", 
                size=0.5,
                alpha=.1) + 
    scale_color_colorblind() +
    scale_x_continuous(breaks=seq(min(mtmelt$day), max(mtmelt$day), 2)) +
    facet_wrap(~symbol, labeller = labeller(symbol=mylabeller),
               nrow=4) + 
    theme_light() + theme(legend.position = "bottom") +
    labs(subtitle=paste(ct,"Old vs Young other top genes"),
         caption=paste('total query (n genes):',dim(g_df)[1] ,
                       '\nGeneSymbol (logFC|padj)'))
  
  # ==================================  plot go ==============
} # end if go_path_doplots
pdf(paste0(resdir,ct,"_GO.pdf"), width=12, height=7)
plot_grid(ggplot() + draw_label(paste("Old vs Young,",ct)),
          plot_grid(gg_go,gg_genenri, ncol=2,labels="AUTO", rel_widths =c(4,5) ) , 
          nrow=2,
          rel_heights = c(1,11) )
gg_genenon
dev.off()
## end
