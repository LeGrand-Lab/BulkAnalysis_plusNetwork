#####  version on Dotplot:
pdf(paste0(resdir,ct,"_genes_OldvsYoung_Dotplot.pdf"))
ggplot(mtmelt, aes(x=day, y=Zscore, group=age , color=age))  +
  geom_point(size=.4, aes(color=age)) + 
  geom_smooth(se=T,
              method="loess", 
              size=0.5,
              alpha=.1) + 
  scale_color_colorblind() +
  facet_wrap(~symbol) + theme_light() +
  labs(title=paste("Old vs Young Top 10 Pathways", ct),
       subtitle="genes with differential dynamics",
       caption=paste("NOTE: no gene was found as significant on padj DESeq interaction
     results \nSelection is based on  *pathway FDR <= 0.05* \n", GG$Description))
dev.off()


####

# pick one gene as example
GG = gp_mod %>% 
  filter(str_detect(Description,"macrophage migration inhibitory")) %>%
  select(Description, FDR, genesEnriched, enriched_terms)

zscogene <- array(zscores[c(GG$genesEnriched), ], 
                  dim= c(length(GG$genesEnriched), dim(zscores)[2]))
colnames(zscogene) <- colnames(zscores)
rownames(zscogene) <- c(GG$genesEnriched)
if(all(rownames(initmeta)==colnames(zscogene))){
  print("rownames(initmeta)==colnames(zscores), ok for reshape2::melt")
  mt = cbind(initmeta,(t(zscogene)))
  mtmelt = reshape2::melt(mt, id.vars=c("newname", "time","type", "age",
                                        "typetimeage"),
                          variable.name="geneid", value.name="Zscore")  
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
  facet_wrap(~geneid) + theme_light()
