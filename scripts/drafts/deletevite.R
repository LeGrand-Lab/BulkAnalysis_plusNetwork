GG = gp_mod %>% 
  filter(!Category=="KEGG" ) %>% 
  slice_min(order_by = FDR, n = 10) %>% 
  select(Description, FDR, genesEnriched, enriched_terms) 

toppathgenes = unique(unlist(str_split(GG$genesEnriched, ","))) #for 




###
#pick genes to plot: genes having padj eq or inf to those found enriched:
candpadj = max(g_df[g_df$id %in% toppathgenes, ]$padj) # max padj enrich genes
viznbtop = dim(g_df %>% filter(padj<=candpadj))[1] + 2 
VIZ <- g_df %>% slice_min(order_by = padj, n = viznbtop)

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
       caption=paste('total query (n genes):',dim(g_df)[1] ,'\n'))

