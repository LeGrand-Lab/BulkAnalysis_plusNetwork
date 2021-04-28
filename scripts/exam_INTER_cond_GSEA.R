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
library(ggrepel) # for labels to points
library("BiocParallel")
register(MulticoreParam(4)) # TODO:  set n of cores depending of available
library(gprofiler2)
library(msigdbr)

setwd("~/BulkAnalysis_plusNetwork/")
source("~/BulkAnalysis_plusNetwork/scripts/gsea_fun.R")
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
dso <- DESeqDataSetFromMatrix(fmat, metadata,
                              design=~ ~ age + time + age:time)

ct = "M1"

all_g_df <- read.table(paste0(resdir, ct, "_INTERagetime.csv"),sep='\t', header = T)

g_df <- all_g_df %>% filter(abs(log2FoldChange)>=1.2)  # this selection for GO
initmeta = metadata %>% filter(type==ct) %>% select(-c(group,sample))
subset <- dso[ unique(g_df$id), dso$type==ct ]
subvsd <- varianceStabilizingTransformation(subset)
zscores = t(scale(t(subvsd@assays@data@listData[[1]]))) 
View(g_df) 
# for M1 there were no significant differences between old young by DESEq2!
# ==================  Gprofiler2 pathways enrichment: ========================

ctgenes <- g_df$padj
names(ctgenes) <- g_df$id
rankedgenes <- sort(ctgenes)
ctres = gprofiler2::gost(names(rankedgenes), organism='mmusculus',
                 ordered_query = TRUE, significant = TRUE,
                 user_threshold = 0.05, evcodes = TRUE,
                 correction_method = "fdr", 
                 sources = c('GO:BP', 'GO:MF', 'GO:CC' , 'KEGG', 
                            'REAC', 'TF', 'MIRNA', 'CORUM', 'HP', 'HPA', 'WP') )
View(ctres$result)
gp_mod = ctres$result %>% select(query, source, term_id,
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
    celltype=ct)
pdf(paste0(resdir,"categALL_",ct,".pdf"), height=10) 
ggplot(gp_mod, 
       aes(celltype,enriched_terms, size = GeneRatio, color=FDR)) +
  geom_point() +
  scale_color_gradient(low="#1AFF1A", high="#4B0092") +
  facet_grid(Category~., scales = "free", space="free") +
  scale_x_discrete(labels = function(x) str_wrap(x, width = 10)) +
  theme(axis.text.y = element_text(size=5)) +
  theme_bw() + 
  labs(title = "Old vs Young, all retreived Categories \n for this cell type",
            caption = paste("Green color most significantly enriched, min FDR:",
                   min(gp_mod$FDR) ) )
dev.off()

## make clearer optimal viz, exclude KEGG
refilter = gp_mod %>%  filter(!Category == "KEGG")
gg_go <- ggplot(refilter, 
       aes( enriched_terms, Count, fill=FDR)) +
  geom_bar(stat="identity", width=0.8) +
  scale_fill_gradient(low="#1AFF1A", high="#4B0092") +
  facet_grid(Category~.,scale="free", space="free") +
  coord_flip() +
  theme(axis.text.y = element_text(size=7)) +
  theme_bw() +
  labs(subtitle = paste("GO/Pathways enrichment"),
       caption = paste("Young vs Old,",ct,
                       "most significant in Greens, min FDR=",min(refilter$FDR),
                       "\nNOTE: explored REACTOME but zero terms retreived" ))


# ################### do for all genes in TOP !! pathways:
GG = gp_mod %>% 
  filter(!Category=="KEGG" ) %>% 
  slice_min(order_by = FDR, n = 10) %>% 
  select(Description, FDR, genesEnriched, enriched_terms) 

toppathgenes = unique(unlist(str_split(GG$genesEnriched, ",")))
zscogene <- array(zscores[c(toppathgenes), ], 
                  dim= c(length(toppathgenes), dim(zscores)[2]))
colnames(zscogene) <- colnames(zscores)
# symbols instead  of ensembl for plotting:
rownames(zscogene) <- genes_df[match(c(toppathgenes),genes_df$Geneid),]$symbol
if(all(rownames(initmeta)==colnames(zscogene))){
  print("rownames(initmeta)==colnames(zscores), ok for reshape2::melt")
  mt = cbind(initmeta,(t(zscogene)))
  mtmelt = reshape2::melt(mt, id.vars=c("newname", "time","type", "age",
                                        "typetimeage"),
                          variable.name="symbol", value.name="Zscore")  
  #mtmelt = full_join(mtmelt, df_clus, by=c("geneid"="id"))
  mtmelt$day = as.numeric(sapply(mtmelt$time, function(x)str_replace(x, "D","")))
  
}else{print("error, zcores colnames and rownames subsetted metadata NOT EQUAL")}

gg_gen <- ggplot(mtmelt, aes(x=day, y=Zscore, group=age , color=age))  +
  geom_point(size=.4, aes(color=age)) + 
  geom_smooth(se=T,
              method="loess", 
              size=0.5,
              alpha=.1) + 
  scale_color_colorblind() +
  facet_wrap(~symbol) + theme_light() +
  labs(subtitle=paste(ct,"Old vs Young representative genes"),
       caption="none of these genes was found  
       \nas significant")


# ==================================  plot go ==============
pdf(paste0(resdir,ct,"_GO.pdf"), width=10, height=6)
plot_grid(ggplot() + draw_label(paste("Old vs Young,",ct)),
plot_grid(gg_go,gg_gen, ncol=2,labels="AUTO") , nrow=2,
rel_heights = c(1,11)
)
dev.off()
# ============================================================================

# =========================== GSEA    ========================================
# ============================================================================
lfcgenes <- all_g_df %>% unique() %>% select(log2FoldChange, id)
lfcgenes$symbols <- genes_df[match(lfcgenes$id,genes_df$Geneid),]$symbol
lfcgenes <- lfcgenes %>% arrange(log2FoldChange)
gseagenes <- lfcgenes$log2FoldChange
names(gseagenes) <- lfcgenes$symbols
gmtneeded = F  # gmt file:
unique(pregmt1$gs_subcat)
if (gmtneeded){
  pregmt1 <- msigdbr(species = "Mus musculus", category = 'C2',
                     subcategory='CP:WIKIPATHWAYS')  
  # [1] "CGP"             "CP:BIOCARTA"     "CP:KEGG"         "CP"           
 # "CP:PID"         "CP:REACTOME"     "CP:WIKIPATHWAYS"
  pregmt2 <- msigdbr(species = "Mus musculus", category = "H")
  thegmt <- bind_rows(pregmt1,pregmt2)
  write.table(thegmt, paste0("stock_gmtfiles/","Hallmark_Wiki.gmt"), sep='\t',
              col.names=T)
}
thegmt <- read.table(paste0("stock_gmtfiles/","Hallmark_Wiki.gmt"), sep='\t',
                     header=T)
#https://davetang.org/muse/2018/01/10/using-fast-preranked-gene-set-enrichment-analysis-fgsea-package/
msigdbr_list = split(x = thegmt$gene_symbol, f = thegmt$gs_name)


fgseaRes <- fgsea::fgsea(pathways = msigdbr_list, 
                      stats = randomlfc,
                      minSize=3,
                      maxSize=Inf, nperm = 100000) 
head(fgseaRes[order(pval), ])

topPathwaysUp <- fgseaRes[][head(order(pval),n=10), pathway]
topPathwaysDown <- fgseaRes[][head(order(pval), n=10), pathway]
topPathways <- c(topPathwaysUp, rev(topPathwaysDown))

uno <- plotGseaTable(msigdbr_list[topPathways], gseagenes, fgseaRes, 
                               gseaParam = 0.5 , render=F)

dos <- plotEnrichment(msigdbr_list[[head(fgseaRes[order(pval), ], 1)$pathway]],
                       gseagenes) + 
              labs(title= ct, subtitle=head(fgseaRes[order(pval), ], 1)$pathway,
       caption="no significant enrichment, M1 genes Old vs Young")

tres <- plotEnrichment(msigdbr_list[[fgseaRes[order(pval), ][7]$pathway]],
                       gseagenes) + 
  labs(title= ct, subtitle=fgseaRes[order(pval), ][7]$pathway, 
       caption="no significant enrichment, M1 genes Old vs Young")

pdf(paste0(resdir,ct,"_GSEA.pdf"), width=10, height=9)
grid.arrange(uno, plot_grid(dos,tres, ncol=2), padding= 20 )
dev.off()


# ============================================================================
##
# https://tex.stackexchange.com/questions/364225/export-tables-from-r-to-latex
