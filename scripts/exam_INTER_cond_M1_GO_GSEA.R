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

ct = "M1"

gologfoldcutoff = 1.2
gopadjcutoff = 1 
topgopath = 10
gsealogfoldcutoff = 0.3
gseapadjcutoff = 1
gseaneeded = F   # set TRUE if needed to run again

prefil_cou <- "data/prefiltered_counts.rds"
metadata.rds <- "data/metadata.rds"
genes_df <- read.table("data/genesinfo.csv",sep="\t",header=T)
# NOTE: color blindness for scales picked from:
#https://davidmathlogic.com/colorblind/#%23D81B60-%231E88E5-%23FFC107-%23004D40

resdir = "exam_INTER_conditions/dynamic/"
fmat <- readRDS(prefil_cou)
metadata <- readRDS(metadata.rds)
metadata <- metadata %>% mutate(typetimeage = paste0(type,".",time,".",age))
##  set which rows to keep :
keep <- apply(fmat, 1, function(row) ifelse(count(row >=5)>= 3, TRUE, FALSE) )
fmat <- fmat[keep,]
# create DESEq2 object to get vst and then zscore 
dso <- DESeqDataSetFromMatrix(fmat, metadata,
                              design= ~ age + time + age:time)

all_g_df <- read.table(paste0(resdir, ct, "_INTERagetime.csv"),sep='\t', header = T)

# =============== prepare data for go and gene viz =================
g_df <- all_g_df %>% filter(abs(log2FoldChange) >= gologfoldcutoff & 
                              padj <= gopadjcutoff ) 
initmeta = metadata %>% filter(type==ct) %>% select(-c(group,sample))
subset <- dso[ unique(g_df$id), dso$type==ct ]
subvsd <- varianceStabilizingTransformation(subset)
zscores = t(scale(t(subvsd@assays@data@listData[[1]]))) 
# ==================================================================
# for M1 there were no significant differences between old young by DESEq2!

# ==================  Gprofiler2 pathways enrichment: ========================
querygenes <- g_df$padj
names(querygenes) <- g_df$id
rankedgenes <- sort(querygenes)
ctres = gprofiler2::gost(names(rankedgenes), organism='mmusculus',
                 ordered_query = TRUE, significant = TRUE,
                 user_threshold = 0.05, evcodes = TRUE,
                 correction_method = "fdr", 
                 sources = c('GO:BP', 'GO:MF', 'GO:CC' , 'KEGG', 
                            'REAC', 'TF', 'MIRNA', 'CORUM', 'HP', 'HPA', 'WP') )
View(ctres$result)
saveRDS(ctres, paste0(resdir,"M1_go_rawobject.rds")) #saved result
## saved result ! last query M1 on gprofiler2  : 16-06-2021 

# prepare for tables and plots:
# -----------------------------------------------------------------
gp_mod = ctres$result %>% filter(precision > 0.01 & recall > 0.01) %>%
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
    enriched_terms = str_trunc(Description, 50, "right"),
    celltype=ct) %>% filter(FDR <= 0.05)

# pdf(paste0(resdir,"categALL_",ct,".pdf")) 
# ggplot(gp_mod, 
#        aes(celltype,enriched_terms, size = GeneRatio, color=FDR)) +
#   geom_point() +
#   scale_color_gradient(low="#1AFF1A", high="#4B0092") +
#   facet_grid(Category~., scales = "free", space="free") +
#   scale_x_discrete(labels = function(x) str_wrap(x, width = 10)) +
#   theme(axis.text.y = element_text(size=5)) +
#   theme_bw() + 
#   labs(title = paste("Old vs Young, all retreived Categories \n",ct),
#       caption = paste(ct, "Old vs Young, pathway enrichment significance (FDR)",
#                 '\nabslfc >=', gologfoldcutoff, 
#                 '  p<=',gopadjcutoff, '  n=',length(querygenes)) )
# dev.off()

## make clearer optimal viz, exclude KEGG (kegg is too much redundance)
gp_mod = gp_mod %>%  filter(!Category == "KEGG")
# 20 relevant pathways by category , to save into table : 
gp_mod= gp_mod %>% group_by(Category)  %>%
  slice_min(order_by = FDR, n = 20) %>% 
  select(Category, ID, Description, FDR, genesEnriched, enriched_terms, 
         query_size, Count, BgRatio, GeneRatio,  celltype ) 
vec_genesymbols <- sapply(gp_mod$genesEnriched, function(x) {
  eids <- unname(unlist(str_split(x, ",")))
  symbs <- genes_df[match(eids, genes_df$Geneid),]$symbol
  return(unname(paste0(symbs, collapse=',')))
})
gp_mod$geneSymbols = vec_genesymbols
write.table(gp_mod %>% select(-enriched_terms), 
            paste0(resdir, "go_gsea_csv/",ct, "_go_filtered.csv"))
# plotting
GG <-  gp_mod %>% group_by(Category) %>%
  slice_min(order_by = FDR, n = topgopath )
gg_go <- ggplot(GG,   aes(enriched_terms, Count, fill=FDR)) +
  geom_bar(stat="identity", width=0.8) +
  scale_fill_continuous(type = "viridis", direction=-1, alpha=.8) +
  facet_grid(Category~.,scale="free", space="free") +
  coord_flip() +
  theme(axis.text.y = element_text(size=7)) +
  theme_bw() +
  labs(subtitle = paste("GO/Pathways enrichment"),
       caption = paste(ct, "Old vs Young, pathway enrichment significance (FDR)",
      '\nabslfc >=', gologfoldcutoff, 
      '  p<=',gopadjcutoff, '  n=',length(querygenes)) )

# do for all genes in TOP !! pathways:
toppathgenes = unique(unlist(str_split(GG$genesEnriched, ","))) 
#pick genes to plot: genes having padj eq or inf to those found enriched:
vizpadj = max(g_df[g_df$id %in% toppathgenes, ]$padj) # max padj enrich genes
viznbtop = dim(g_df %>% filter(padj<=vizpadj))[1] + 2 
VIZ <- g_df %>% slice_min(order_by = padj, n = viznbtop)
VIZ$symbol <- genes_df[match(VIZ$id,genes_df$Geneid),]$symbol
VIZ$wasenriched <- ifelse(VIZ$id %in% toppathgenes, 1, 0)

# build labeller
mylabeller <- paste0(VIZ$symbol," ", " (",
                          round(VIZ$log2FoldChange,2),"|",
                     round(VIZ$padj,2),")" )
names(mylabeller) <- VIZ$symbol

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
  facet_wrap(~symbol, labeller = labeller(symbol=mylabeller)) + 
  theme_light() + theme(legend.position = "bottom") +
  labs(subtitle=paste(ct,"Old vs Young, GO/pathways retrieved"),
       caption=paste('total query (n genes):',dim(g_df)[1] ,
                     '\nGeneSymbol (logFC|padj)'))

gg_genenon <- ggplot(mtmelt %>% filter(wasenriched==0), 
                    aes(x=day, y=Zscore, group=age , color=age))  +
  geom_point(size=.4, aes(color=age)) + 
  geom_smooth(se=T,
              method="loess", 
              size=0.5,
              alpha=.1) + 
  scale_color_colorblind() +
  facet_wrap(~symbol, labeller = labeller(symbol=mylabeller)) + theme_light() +
  labs(subtitle=paste(ct,"Old vs Young representative genes"),
       caption=paste('total query (n genes):',dim(g_df)[1] ,
                     '\nGeneSymbol (logFC|padj)'))

# ==================================  plot go ==============
pdf(paste0(resdir,ct,"_GO.pdf"), width=10, height=6)
plot_grid(ggplot() + draw_label(paste("Old vs Young,",ct)),
plot_grid(gg_go,gg_genenri, ncol=2,labels="AUTO") , nrow=2,
rel_heights = c(1,11) )
gg_genenon
dev.off()



# =========================== GSEA    ========================================
# https://davetang.org/muse/2018/01/10/using-fast-preranked-gene-set-enrichment-analysis-fgsea-package/
# https://bioconductor.org/packages/release/bioc/vignettes/fgsea/inst/doc/fgsea-tutorial.html
# https://rdrr.io/bioc/fgsea/f/vignettes/fgsea-tutorial.Rmd
# ============================================================================
gmtneeded = F  # gmt file:
if (gmtneeded){
  pregmt1 <- msigdbr(species = "Mus musculus", category = 'C2',
                     subcategory=c('CP:REACTOME')  )
  # [1] "CGP"             "CP:BIOCARTA"     "CP:KEGG"         "CP"           
  # "CP:PID"         "CP:REACTOME"     "CP:WIKIPATHWAYS"
  pregmt2 <- msigdbr(species = "Mus musculus", category = "H") #HALLMARK
  thegmt <- bind_rows(pregmt1,pregmt2)
  write.table(thegmt, paste0("stock_gmtfiles/","Hallmark_React.gmt"), sep='\t',
              col.names=T)
}

if(gseaneeded){
  lfcgenes <- all_g_df %>% filter(abs(log2FoldChange) >= gsealogfoldcutoff &
                                    padj <= gseapadjcutoff) %>% 
    arrange(padj) %>% distinct(id, .keep_all=T) %>%
    select(log2FoldChange, id) 
  
  lfcgenes$symbols <- genes_df[match(lfcgenes$id,genes_df$Geneid),]$symbol
  lfcgenes <- lfcgenes %>% arrange(log2FoldChange)
  gseagenes <- lfcgenes$log2FoldChange
  names(gseagenes) <- lfcgenes$symbols
  barplot(sort(gseagenes,decreasing=T))
  
  thegmt <- read.table(paste0("stock_gmtfiles/","Hallmark_React.gmt"), sep='\t',
                       header=T)
  
  msigdbr_list = split(x = thegmt$gene_symbol, f = thegmt$gs_name)
  set.seed(42)
  fgseaRes <- fgsea::fgsea(pathways = msigdbr_list, 
                           stats = gseagenes,
                           minSize=3,
                           maxSize=Inf, nperm = 100000) 
  head(fgseaRes[order(pval), ])
  saveRDS(fgseaRes, file=paste0(resdir,ct,'_GSEAdatafull.rds'))
  
  topPathwaysUp <- fgseaRes[ES>0][head(order(padj),n=10), pathway]
  topPathwaysDown <- fgseaRes[ES<0][head(order(padj), n=10), pathway]
  topPathways <- c(topPathwaysUp, rev(topPathwaysDown))
  uno <- plotGseaTable(msigdbr_list[topPathways], gseagenes, fgseaRes, 
                       gseaParam = 0.5 , render=F)
  dos <- plotEnrichment(msigdbr_list[[topPathwaysUp[1]]],    gseagenes) + 
    labs(title= ct, subtitle=topPathwaysUp[1],
         caption=paste(ct,"Old vs Young, abslfc >=",gsealogfoldcutoff,
                       "p<=",gseapadjcutoff, "n=", length(gseagenes) ))
  tres <- plotEnrichment(msigdbr_list[[topPathwaysDown[1]]], gseagenes) + 
    labs(title= ct, subtitle=topPathwaysDown[1],
         caption=paste(ct,"Old vs Young, abslfc >=",gsealogfoldcutoff,
                       "p<=",gseapadjcutoff, "n=", length(gseagenes) ))
  pdf(paste0(resdir,ct,"_GSEA.pdf"), width=10, height=9)
  grid.arrange(uno, plot_grid(dos,tres, ncol=2), padding= 20 )
  dev.off()
  
  # save table containing top listed pathways and matching genes :
  # ----------------------------------------------------------------------
  # fgseaRes <- readRDS(paste0(resdir,ct,'_GSEAdatafull.rds'))
  fgsea_top <- fgseaRes[padj <= 0.3]
  fgsea_top <- fgsea_top %>% mutate(listOfgenes = 
                        str_replace_all(leadingEdge, c('c\\("' = '', 
                                              '"' = '',
                                              '\\)' = ''))) %>% 
    relocate(listOfgenes, .after = size) %>%
    mutate(regulation = case_when(ES > 0 ~ "Up", TRUE ~ "Down")) %>%
    relocate(regulation, .after = ES) %>% select(-c(nMoreExtreme,leadingEdge,size))
  
  nbofgenes = sapply(fgsea_top$listOfgenes, 
                     function (x){ 
                       elems = unname(str_split(x, ', ')[[1]])
                       length(unname(elems))  })
  fgsea_top$ngenes = nbofgenes
  
  write.table(fgsea_top, paste0(resdir, "go_gsea_csv/", ct, '_GSEA_topPaths.csv'),
              sep="\t", col.names = T, row.names = F)
  
  # ----------------------------------------------------------------------
  # end saving table
  collapseve = F # this 'collapse' was very useful, check ..collapsed.csv
  if (collapseve = T){
    collapsedPathways <- collapsePathways(fgseaRes[order(pval)][padj < 0.01], 
                                          msigdbr_list, gseagenes)
    mainPathways <- fgseaRes[pathway %in% collapsedPathways$mainPathways][
      order(-NES), pathway]
    plotGseaTable(msigdbr_list[mainPathways], gseagenes, fgseaRes, 
                  gseaParam = 0.5)
    # saving table : 
    mainPathTab <- fgseaRes[pathway %in% collapsedPathways$mainPathways][order(-NES),]
    tosave =  tibble("num"=names(unlist(mainPathTab[,1])))
    for (i in names(mainPathTab)){
      print(i)
      tmp = unlist(mainPathTab[[i]])
      if (i == 'leadingEdge'){
        tmp = paste0(unlist(mainPathTab[[i]]), collapse = ",")
      }
      tosave <- cbind(tosave, unname(tmp))
    }
    tosave$num <- NULL
    colnames(tosave) <- names(mainPathTab)
    write.table(tosave,paste0(resdir,"notrust/",ct,"_GSEAcollapsed.csv"), sep='\t',
              col.names = T, row.names = T )
  } #end if collapseve
}

# ============================================================================
##
# https://tex.stackexchange.com/questions/364225/export-tables-from-r-to-latex
