# https://stackoverflow.com/questions/52970275/ggalluvial-order-flow-of-lines-based-on-a-variable-within-stratum
library(ggalluvial)
data(majors)
majors$curriculum <- as.factor(majors$curriculum)
ggplot(majors,
       aes(x = semester, stratum = curriculum, alluvium = student,
           fill = curriculum, label = curriculum)) +
  scale_fill_brewer(type = "qual", palette = "Set2") +
  geom_flow(stat = "alluvium", lode.guidance = "frontback",
            color = "darkgray") +
  geom_stratum() +
  theme(legend.position = "bottom") +
  ggtitle("student curricula across several semesters")
malluv <- to_alluvia_form(majors,
                          key="semester", value="curriculum",
                          id = "student")
gg <- ggplot(malluv,
             aes(axis1 = CURR1, axis2 = CURR7, axis3 = CURR13))
gg +
  geom_alluvium(aes(fill = as.factor(student)), width = 2/5, discern = TRUE) +
  geom_stratum(width = 2/5, discern = TRUE) +
  geom_text(stat = "stratum", discern = TRUE, aes(label = after_stat(stratum)))
#=========================================

NEWMELTED <- givemeltedDataframe(pathsFiltered, 2, 15, "REACTOME")
# print("adding information about genes (padj, lfc), from 'rds/shot_rds_full.rds'")
# toto <- NEWMELTED %>% select(unigene, gene, celltype, day)
# fuu <- fullDEsta %>% select(symbol, padj, log2FoldChange, type, day, baseMean) %>% 
#   mutate(gene = symbol, celltype = type, genepadj = padj)
# tototo <- left_join(toto, fuu, by = c("day", "gene", "celltype" ) )
# tototo <- tototo  %>% select(unigene, genepadj, log2FoldChange, day, baseMean )
# ttp = left_join(NEWMELTED, tototo, by=c("unigene", "day"))  %>% unique() 
# NEWMELTED = ttp
# rm(toto, fuu, tototo, ttp)

gogo = NEWMELTED %>% group_by(pathway, day, celltype) %>% 
  summarize(pathway,day, celltype, ngenes=n(), padj, NES) %>% unique()
# my way to select some paths
mysele = gogo %>% group_by(pathway) %>% 
  summarize(pathway, noccu=n(), padj) %>% 
  unique() %>% filter(noccu >=5 , padj <= 0.005) %>% pull(pathway)

#allcases = NEWMELTED %>% group_by(pathway, day, celltype) %>% unique() %>%
 # filter(pathway %in% mysele)

## ::: ! 

dummy <- data.frame(gogo %>% filter(pathway %in% mysele) %>%
                      mutate(pathway = sapply( pathway, function(x) str_replace(x, "REACTOME_",""))))
dummy$pathway <- as.factor(dummy$pathway)

dummy <- dummy %>% filter(pathway %in% c("CELL_CYCLE_CHECKPOINTS", "CELLULAR_SENESCENCE"))

View(dummy)
dummy <- rbind(dummy, c("CELL_CYCLE_CHECKPOINTS", "D0", "none",1,1, -0.1))
dummy <- rbind(dummy, c("CELLULAR_SENESCENCE", "D2", "M1", 1, 1, -0.1 ))
dummy <- rbind(dummy, c("CELLULAR_SENESCENCE", "D7", "M1", 1, 1, -0.1 ))
dummy <- rbind(dummy, c("CELLULAR_SENESCENCE", "D7", "ECs", 1,1, -0.1))

dummy2 <- dummy %>% filter(day %in% c("D4","D7") & celltype %in% c("M1","ECs") &
                             pathway == "CELLULAR_SENESCENCE")
dummy1 <- dummy %>% filter(day %in% c("D4","D7") & celltype %in% c("M1","ECs") &
                             pathway == "CELL_CYCLE_CHECKPOINTS")
dummy3 <- dummy %>% filter(day %in% c("D2", "D4","D7") & 
                             pathway == "CELLULAR_SENESCENCE")

fifo = rbind(dummy1,dummy2) %>% select(pathway,day, celltype)
ggplot(fifo, 
       aes(x = day, 
           stratum = pathway,
           alluvium = celltype,
           fill = celltype,
           label = pathway)) + 
  scale_fill_manual(labels = c("ECs","M1", "M2", "sCs", "none"), 
                    values = c("blue", "red", "violet","royalblue", "whitesmoke")) +  
  geom_flow(stat = "alluvium", "lode.guidance" = "frontback",
            color = "darkgray") + 
  geom_stratum(alpha=0.3) +
  geom_text(stat = "stratum", aes(label=after_stat(stratum)))+
  theme(legend.position="bottom") + 
  ggtitle("Old vs Young differentially enriched pathways (GSEA")


ggplot(dummy2, 
       aes(x = day, y = NES,
           stratum = pathway,
           alluvium = celltype,
           fill = celltype,
           label = pathway)) + 
  scale_fill_manual(labels = c("ECs","M1", "M2", "sCs", "none"), 
                    values = c("blue", "red", "violet","royalblue", "whitesmoke")) +  
  geom_flow() + 
  geom_stratum(alpha=0.3) +
  geom_text(stat = "stratum", aes(label=after_stat(stratum)))+
  theme(legend.position="bottom") + 
  ggtitle("Old vs Young differentially enriched pathways (GSEA")



#=========
  
ss <- list()
for (d in c('D0', 'D2', 'D4', 'D7')){  
  moo = sapply(colnames(dummy), function(x) paste0(x,d))
  ss[[d]] <- dummy %>% filter(day == d)
  colnames(ss[[d]]) <- moo
}
ouf <- full_join(ss[['D2']], ss[['D4']], by=c("pathwayD2"="pathwayD4", "celltypeD2"="celltypeD4"))
ggplot(ouf,
       aes(axis1=pathwayD2, axis2=pathwayD2, y = NES)) +
  geom_alluvium(aes(fill=celltype))
         
#=========
#=========
#=========

tiy <- sample_n(xxxxxxxxxxxxxx, 10)

CUTOFFPATHWAY = 0.8

vv <- NEWMELTED %>% filter(padj <= CUTOFFPATHWAY) %>%
       mutate(pathway = sapply( pathway, function(x) str_replace(x, "REACTOME_","")))
table(vv$celltype)

gogo = NEWMELTED %>% group_by(pathway, day, celltype) %>% 
  summarize(pathway,day, celltype, ngenes=n(), padj, NES) %>% unique()

ggplot(data = gogo, 
       aes(x=pathway, y=-log10(padj))) + 
  geom_point(size = 3, colour = "royalblue") +
  geom_segment(aes(xend=pathway, yend = 0), size= 0.5) + 
  theme(axis.text.x = element_text(size = 2, angle = 90)) + coord_flip()

library(viridis)
library(cowplot)

bycell <- ggplot(data = gogo ) + 
  geom_jitter( aes(x=-(padj), y=NES, color = celltype, size=ngenes) , alpha = 0.7) +
  scale_color_manual(values = c("darkgray", "gold", "red", "violet", "green", "royalblue")) +
  theme_light()
byday <- ggplot(data = gogo ) + 
  geom_jitter( aes(x=-(padj), y=NES, color = day) , alpha = 0.7) +
  scale_color_viridis_d(option="rocket", direction=-1) +
  theme_light()

bycell2 <- ggplot(data = gogo) + 
  geom_jitter( aes(x=-log10(padj), y=NES, color = celltype, size=ngenes) , alpha = 0.7) +
  scale_color_manual(values = c("darkgray", "gold", "red", "violet", "green", "royalblue")) +
  theme_light( )
byday2 <- ggplot(data = gogo) + 
  geom_jitter( aes(x=-log10(padj), y=NES, color = day) , alpha = 0.7) +
  scale_color_viridis_d(option="rocket", direction=-1) +
  theme_light()


pdf(paste0(odir, "GSEA/GSEAbyctbyday_panorama.pdf"), width = 12, height=10)
title <- ggdraw() + draw_label("Pathways from GSEA result", size = 14)
caption <- ggdraw() +
  draw_label("GSEA was performed separately by celltype and by day\n No Tau filtering",
             size = 9)
plot_grid(title, plot_grid(bycell, byday, bycell2, byday2), caption, ncol= 1,
          rel_heights = c(0.2,1, 0.1))  
dev.off()
`



     
