# GO , REACTOME, KEGG and other functional databases enrichment
# requires results obtained from 'dynamics_intra_Spec'
# Enrichment on tissue-Specific-highly-dynamic genes only

library(dplyr)
library(tidyverse)
require(gridExtra) # calling grid.arrange
library("BiocParallel")
register(MulticoreParam(4)) # TODO:  set n of cores depending of available
library(DESeq2)
library(cluster)
library(reshape2)
library(RColorBrewer)
library(ggthemes)
library(cowplot)
library(gprofiler2)

setwd("~/BulkAnalysis_plusNetwork/")
enrdir = "dynamicsIntra_Spec/Gprofiler2_res/"
here.age = "Young"  # our "reference"
print("setting top n (lowest p.adjust) terms to select")
nbtop <- 10
print("setting databases source")
pickcat <- c("GO:BP",  "GO:MF" , "REAC" )
cellcolors = list(
  "ECs"="#0072B2",
  "FAPs"="#F0E442",
  "M1" = "#D55E00",
  "M2" =  "#CC79A7",
  "Neutro" =  "#009E73",
  "sCs" = "#56B4E9" 
)
topTerms <- list()
topPlots <- list()
#sign_list = readRDS(paste0("dynamicsIntra_Spec/",here.age, "IntraDyn_list.rds"))
# exclude "Neutro" to avoid errors
for(cty in names(cellcolors)[!names(cellcolors)=="Neutro"]){
  clus_df = read.table(paste0("dynamicsIntra_Spec/",
             "clustersKmeans_",here.age,"_",cty,".csv"), sep='\t', header=T)
  i.clusts = sort(unique(clus_df$cluster))
  # make a list of clusters, for GO analysis
  i.cllist = list()
  for(k in i.clusts){
    print(k)
    i.cllist[[k]] = clus_df %>% filter(cluster==k) %>% pull(id)
  }
  i.res = gost(query = i.cllist, organism="mmusculus")
  gp_mod = i.res$result %>% select(query, source, term_id,
                              term_name, p_value,query_size, 
                              intersection_size,term_size,
                    effective_domain_size)
  print("if gene ids wanted, check option 'evcodes', but will run slower")
  gp_mod$GeneRatio = gp_mod$intersection_size/gp_mod$query_size
  gp_mod$BgRatio = gp_mod$term_size/gp_mod$effective_domain_size
  names(gp_mod) = c("Cluster", "Category", "ID", "Description", "FDR",
   "query_size", "Count", "term_size", "effective_domain_size",
    "GeneRatio", "BgRatio")
  write.table(gp_mod, paste0(enrdir, here.age, "_",cty,"_gprofiler2_res.csv"),
              col.names=T, row.names = F)
  print("these are the terms/databases we can pick results from")
  print(unique(gp_mod$Category))  #   # Visualize across ALL enriched terms/databases, testing relevance:
  relevantEnr = gp_mod %>% filter(Category %in% pickcat & FDR <= 0.05 ) %>%
    group_by(Cluster, Category) %>% 
    slice_min(order_by = FDR, n = 5) %>% 
    mutate(enriched_terms=str_trunc(Description, 50, "right"))
}

allTOP.here.age <- bind_rows(topTerms)
write.table(allTOP.here.age %>% mutate(age=here.age),
            file=paste0(enrdir, here.age, "TOPtermsSummary.csv"),
            col.names = T, row.names = F)
for(cty in names(cellcolors)[!names(cellcolors)=="Neutro"]){
  relevantEnr <- allTOP.here.age %>% filter(Celltype==cty)
  topPlots[[cty]] <- ggplot(relevantEnr, aes(Cluster, enriched_terms, 
                                             size = GeneRatio, color=FDR)) +
    geom_point() +
    scale_color_gradient(low=cellcolors[[cty]], high="gray11") +
    theme(axis.text.y=element_text(size=6)) + 
    facet_grid(Category~., scale = "free_y") +
    labs(title=paste(cty, here.age, "tissue-specific-dynamic-genes"),
         caption=paste("terms filtered by FDR <= 0.05, and top",nbtop)) +
    theme_light() 
  topTerms[[cty]] <- relevantEnr %>% mutate(CellType=cty)
  
}

pdf(paste0(enrdir, here.age,"_", names(topPlots)[1] ,"_plot.pdf"), width=8, height=7)
topPlots[1]
dev.off()

pdf(paste0(enrdir, here.age, "_", names(topPlots)[2],"_plot.pdf"), width=8, height=7)
topPlots[2]
dev.off()

pdf(paste0(enrdir, here.age, "_",names(topPlots)[3],"_plot.pdf"), width=6, height=4)
topPlots[3]
dev.off()

pdf(paste0(enrdir, here.age, "_",names(topPlots)[4],"_plot.pdf"), width=7,  height=4)
topPlots[4]
dev.off()

pdf(paste0(enrdir, here.age, "_",names(topPlots)[5],"_plot.pdf"), width=8, height=7)
topPlots[5]
dev.off()

## ======================Old run: exploratory across all databases results
# Visualize across ALL enriched terms/databases, testing relevance:
testneeded = F  # keep FALSE as already performed, no longer needed
if (testneeded){
  cty="FAPs"  # or whatever
  testtopn = 5
  clus_df = read.table(paste0("dynamicsIntra_Spec/",
                              "clustersKmeans_",here.age,"_",cty,".csv"), sep='\t', header=T)
  i.clusts = sort(unique(clus_df$cluster))
  # make a list of clusters, for GO analysis
  i.cllist = list()
  for(k in i.clusts){
    print(k)
    i.cllist[[k]] = clus_df %>% filter(cluster==k) %>% pull(id)
  }
  i.res = gost(query = i.cllist, organism="mmusculus")
  gp_mod = i.res$result %>% select(query, source, term_id,
                                   term_name, p_value,query_size, 
                                   intersection_size,term_size,
                                   effective_domain_size)
  gp_mod$GeneRatio = gp_mod$intersection_size/gp_mod$query_size
  gp_mod$BgRatio = gp_mod$term_size/gp_mod$effective_domain_size
  names(gp_mod) = c("Cluster", "Category", "ID", "Description", "FDR",
                    "query_size", "Count", "term_size", "effective_domain_size",
                    "GeneRatio", "BgRatio")
  print("all gprofiler2 consulted databases: ")
  print(unique(gp_mod$Category))
  #[1] "GO:MF" "KEGG"  "REAC"  "WP"    "HP"    "MIRNA" "GO:CC" "TF"   
  test = gp_mod %>% group_by(Cluster, Category) %>% 
    slice_min(order_by = FDR, n = testtopn)
  pdf(paste0(enrdir,"windowsALL4Test_",cty,".pdf"), width=15, height=18)
  ggplot(test, aes(Cluster,Description, size = GeneRatio, color=FDR)) +
    geom_point() +
    scale_color_continuous() +
    facet_grid(Category~.)  +
    scale_x_discrete(labels = function(x) str_wrap(x, width = 10)) +
    theme(axis.text.y = element_text(size=.6))
  theme_linedraw()
  dev.off()
}# end if test
  

