ct = "M1"
mat.ct <- fmat[,str_detect(colnames(fmat), ct)]
meta.ct <- metadata %>% filter(str_detect(rownames(metadata),ct))
x.keep <- apply(mat.ct, 1, function(row) ifelse(sum(row >=5)>= 3, TRUE, FALSE))
mat.ct <- mat.ct[x.keep,]
ds.o <- DESeqDataSetFromMatrix(mat.ct, meta.ct, design = ~ age + time + age:time)
ds.o$age <- relevel(ds.o$age, ref="Young")
timeordered <- sort(unique(ds.o$time))
ds.o$time <- factor(ds.o$time, levels=timeordered)
d.r <- DESeq(ds.o, full=~age + time + age:time)
print("ended DE step, picking interaction contrast")
# the interaction term, answering: is the time effect *different* across ages ?
mycontrast <- paste0("ageOld.time", timeordered[2])
tint <- results(d.r, name=mycontrast)
print(head(tint))
tint$id <- rownames(tint)
tint <- as_tibble(tint) %>% 
  mutate(contrast=paste0(rev(timeordered),collapse="_vs_"),
         celltype = ct ) 
write.table(tint %>% filter(abs(log2FoldChange)>= 1.2 &
                              padj <= 0.1), # see non significant also
            file=paste0(resdir, ct, "V2_INTERagetime.csv"),
            sep='\t', col.names=T, row.names = F)
print(paste("saved csv file for INTERaction, ",ct))




# ====== for spaghetti plot ================

filteredDF <- read.table(paste0(resdir,"filteredDynamicLRT.csv"),
                         sep='\t',header=T)
entire.ds <- DESeqDataSetFromMatrix(countData = fmat,
                                    colData = metadata,
                                    design = ~ age + time + age:time )
prep = list()
for(ct in names(allresu_l)){
  print(ct)
  selectgenes <- unique(filteredDF[filteredDF$type==ct, ]$id)
  # get rlog counts for this genes for this celltype
  chosen <- entire.ds[rownames(entire.ds) %in% selectgenes, entire.ds$type==ct]
  chosen <- rlog(chosen)
  tmp <- data.frame("id"=rownames(chosen)) # tmp dataframe to stock stuff
  disptimes <- sort(unique(chosen$time))
  for (thisage in c("Young","Old")){
    for (thistime in disptimes){
      tmptmp <- chosen[,chosen$age==thisage & chosen$time==thistime]
      vecmedian <- apply(assay(tmptmp),1, function(row) median(row))
      vecmean <- apply(assay(tmptmp),1, function(row) mean(row))
      tmp[[paste0(thisage,".",ct,".",thistime,".median")]] <- vecmedian
      tmp[[paste0(thisage,".",ct,".",thistime,".mean")]] <- vecmean
    }
  }
  prep[[ct]] <- tmp
}#end  
saveRDS(prep,file=paste0(resdir,"dyn_prep_spaghetti.rds") ) 


aa_l <- list()
# take only medians for this plot:
for(ct in names(prep)){
  tmp <- prep[[ct]][,str_detect(colnames(prep[[ct]]),".median")] # take only medians
  chosengenes <- prep[[ct]]$id
  chosencols <- colnames(tmp)
  dfdf <- as.data.frame(t(as.matrix(tmp))) #transpose
  colnames(dfdf) <- chosengenes
  dfdf$detailcond <- rownames(dfdf) 
  dfdf$age <- sapply(dfdf$detailcond, function (x) str_split(x,"\\.")[[1]][1])
  dfdf$celltype <- sapply(dfdf$detailcond, function (x) str_split(x,"\\.")[[1]][2])
  dfdf$time <- sapply(dfdf$detailcond, function (x) str_split(x,"\\.")[[1]][3])
  aa <- melt(dfdf,time=time, age=age, celltype=celltype, median_rlogcount=value)
  aa$age <- factor(aa$age)  
  aa$median_rlogcount <- aa$value 
  aa <- aa %>% mutate(time_day=case_when(
    time == "D0" ~ 0,
    time == "D2" ~ 2,
    time == "D4" ~ 4,
    time == "D7" ~ 7))
  vecsignifi = allresu_l[[ct]][allresu_l[[ct]]$padj <= 0.05,]$id                                                                         
  aa <- aa %>% mutate(significant=ifelse(variable %in% vecsignifi,
                                         "FDR <= 0.05","No Sig"))
  aa_l[[ct]] <- aa  
}


bigaa <- dplyr::bind_rows(aa_l)
col_vir <- viridis_pal(begin=0,end=1)(10)  # nice scale, to pick from

pdf(paste0(resdir,"dyn_consecutivetp.pdf"), width=14, height = 6)  
ggplot(data = bigaa %>% mutate(gene_age=paste0(variable,age)), 
       aes(x=as.numeric(time_day), y=median_rlogcount, 
           group = gene_age,
           color=age)) + 
  geom_point( aes(color=age), alpha=.3, size=.5) + 
  geom_line( aes(linetype=significant, color=age),alpha=.3, size=.8) +
  scale_color_manual(values=c(col_vir[2], col_vir[6])) +
  scale_y_log10() + labs(xlab="time_point (day)") +
  theme_bw() + facet_grid(~celltype) + labs(title="Dynamics of expression",
                                            x="time_point (day)",
                                            caption="only genes which after LRT result in FDR < 0.1 
            for at least one consecutive time point comparison")
dev.off()
## plot only significant:
bigaaSIG <- bigaa %>% filter(significant=="FDR <= 0.05")
pdf(paste0(resdir,"dyn_consecutivetpSIGNIF.pdf"), width=14, height = 6)  
ggplot(data = bigaaSIG %>% mutate(gene_age=paste0(variable,age)), 
       aes(x=as.numeric(time_day), y=median_rlogcount, 
           group = gene_age,
           color=age)) + 
  geom_point( aes(color=age), alpha=.3, size=.5) + 
  geom_line( aes( color=age),alpha=.3, size=.8) +
  scale_color_manual(values=c(col_vir[2], col_vir[6])) +
  scale_y_log10() + 
  theme_bw() + facet_grid(~celltype) + 
  labs(title="only significant Genes (by FDR <= 0.05)",
       x="time_point (day)")
dev.off()








