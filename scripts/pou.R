#========= this draft requires running dynamics_intra.R and keep Environment !!
# joha GL
#
print("this draft requires running dynamics_intra.R and keep Environment !!")
print(paste(
   "(!)pdf files 'ECs_arbitrary9_geompath' and 'ECs_arbitrary9_smooth' were",
   "produced using arbitrary 9 clusters, the plotted data is the same in both.",
   "Some clusters are very similar!!, better use silhouette criteria!")
)

opticlust <- list( "ECs" =9) # arbitrary

here.age = "Old"
cty <- "ECs"
print(paste("dynamicAndSpecific genes found in **", here.age, "** mice"))
DYNSPE = readRDS(paste0(resdir, here.age,"IntraDyn_list.rds"))
for(i in names(DYNSPE)){
   print(paste(i, "==> ", length(unique(DYNSPE[[i]]$id))))
}
# [1] "M1 ==>  148"
# [1] "M2 ==>  261"
# [1] "FAPs ==>  911"
# [1] "sCs ==>  992"
# [1] "ECs ==>  476"

agemat = fmat[,str_detect(colnames(fmat), here.age)]
keep <- apply(agemat, 1, function(row) ifelse(sum(row >=5)>= 3, TRUE, FALSE) )
agemat <- agemat[keep,]
here.meta = metadata %>% filter(age==here.age)
dso <- DESeqDataSetFromMatrix(agemat, here.meta, design=~TYPEtime)
subset <- dso[ unique(DYNSPE[[cty]]$id), dso$type==cty ]
subvsd <- varianceStabilizingTransformation(subset)
zscores = t(scale(t(subvsd@assays@data@listData[[1]]))) 

kclushere = kmeans(zscores, centers=opticlust[[cty]])
df_clus = tibble("id"=names(kclushere$cluster), "cluster"=kclushere$cluster)
colors = c("pink","green","lightblue","gold","gray", "blue","cyan","cadetblue","purple")
p <- ggplot()
aplist <- list(p, p,p,p,p,p,p,p,p)
initmeta <- here.meta %>% filter(type==cty)
## with geom_path
for (i in 1:475){
   fly = initmeta
   fly$Z_score = zscores[rownames(zscores)[i],]
   cluster = df_clus$cluster[i]
   fly$day = as.numeric(sapply(fly$time, function(x)str_replace(x, "D","")))
   print(fly)
   aplist[[cluster]] <- aplist[[cluster]] + geom_path(data=fly, 
            mapping=aes(x =day, y=Z_score), col=colors[cluster])
}
plot_grid(plotlist=aplist, nrow=2)
## with geom_smooth
yplo = list(p, p,p,p,p,p,p,p,p)
for (i in 1:475){
   fly = initmeta
   fly$Z_score = zscores[rownames(zscores)[i],]
   cluster = df_clus$cluster[i]
   fly$day = as.numeric(sapply(fly$time, function(x)str_replace(x, "D","")))
   print(fly)
   yplo[[cluster]] <- yplo[[cluster]] + 
      geom_smooth(data=fly, mapping=aes(x =day, y=Z_score), 
                  col=colors[cluster], alpha=.2)
}
plot_grid(plotlist=yplo, nrow=2)

here.age = "Young"
pp = readRDS(paste0(resdir,"complexGGplot_trends",here.age,".rds"))
pdf(paste0(resdir,"SpecificDyn_",here.age,".pdf"), width=7, height=14)
plot_grid(plotlist = pp, ncol= 1, rel_heights = c(2,2,1,2,2))
dev.off()
