## 
# Dynamics of diff expression (in time)
# Each cell type separately (Neutro EXCLUDED as single-time-point measured)
# --
# Joha GL
##
#library(reticulate)
#use_condaenv(condaenv="~/prog_bio/anaconda3/envs/MYVENV",
#             conda="~/prog_bio/anaconda3/bin/conda")
library(dplyr)
library(tidyverse)
library(ggplot2)
library(cowplot)
library(forcats)
library(RColorBrewer)
library(MASS)
library(pheatmap)
library(DESeq2)
library(ggsci) # publishing palettes
library(cowplot)
library(gridExtra)
library(reshape2)
library(ggthemes)
library(scales) # 'viridis_pal' et autres pal
library(ggrepel) # for labels to points
library("apeglm") # BiocManager::install("apeglm")
library("BiocParallel")
library("ashr") # install.packages("ashr")
register(MulticoreParam(4)) # TODO:  set n of cores depending of available
library(pheatmap)


setwd("~/bulk_analysis/")
prefil_cou <- "data/prefiltered_counts.rds"
metadata.rds <- "data/metadata.rds"
genes_df <- read.table("data/genesinfo.csv",sep="\t",header=T)

resdir = "plotsDE/"
DEdir = "signaturestypes/"
fmat <- readRDS(prefil_cou)
metadata <- readRDS(metadata.rds)
metadata <- metadata %>% mutate(timetype=paste0(time,".",type)) 

keep <- apply(fmat, 1, function(row) ifelse(count(row >=5)>= 3, TRUE, FALSE) )
fmat <- fmat[keep,]

allresu_l <- list() #to stock all results
allct <- sort(unique(metadata$type))
# loop all celltypes (excluding Neutro as it has only one timepoint)
allct <- allct[!allct=="Neutro"]
for (ct in allct){  
  mat.ct <- fmat[,str_detect(colnames(fmat), ct)]
  meta.ct <- metadata %>% filter(str_detect(rownames(metadata),ct))
  x.keep <- apply(mat.ct, 1, function(row) ifelse(sum(row >=5)>= 3, TRUE, FALSE))
  mat.ct <- mat.ct[x.keep,]
  ds.o <- DESeqDataSetFromMatrix(countData = mat.ct,
                                 colData = meta.ct,
                                 design = ~ age + time + age:time )
  colData(ds.o)$age <- relevel(colData(ds.o)$age, ref="Young")
  # detect consecutive points for current dataset
  allpoints <- as.character(levels(colData(ds.o)$time))
  tmp_res_consec = list()
  for (k in 1:(length((allpoints))-1)){
    consec <- c(allpoints[k],allpoints[k+1])
    print(consec)
    ds.consec <- ds.o[,as.character(ds.o$time) %in% consec]
    colData(ds.consec)$time <- factor(colData(ds.consec)$time,
                                      levels=sort(unique(colData(ds.consec)$time)))
    ddsTC <- DESeq(ds.consec, test="LRT", reduced = ~ age + time)
    # justification: see all this https://support.bioconductor.org/p/114179/
    resultsNames(ddsTC)
    resTC <- results(ddsTC, contrast = c("age", "Old", "Young"), parallel=T)
    resTC$id <- rownames(resTC)
    # select no stringent FDR , to keep most genes as possible in this step
    tmp_res_consec[[k]]  <- as_tibble(resTC) %>% 
      arrange(padj)  %>%
      filter(padj < 0.1) %>% mutate(type=ct, 
                                    contrast = paste0(consec, collapse="vs"))
  } # end for each consecutive points
  allresu_l[[ct]] <- dplyr::bind_rows(tmp_res_consec) 
}# end for each cell type

# make sure all results successfully stocked into list (only view)
View(as_tibble(allresu_l[["ECs"]]) %>% filter(padj < 0.05))
View(as_tibble(allresu_l[["FAPs"]]) %>% filter(padj < 0.05))
View(as_tibble(allresu_l[["M1"]]) %>% filter(padj < 0.05))
View(as_tibble(allresu_l[["M2"]]) %>% filter(padj < 0.05))
View(as_tibble(allresu_l[["sCs"]]) %>% filter(padj < 0.05))

# significant are too few  for M1 (2) and M2 (28)
# NOTE: there are genes in common between two time contrasts for given cell type
#       but most dyn exp are captured only at one time contrast
# so we still keep all, add gene symbols ??? and save into rds file
saveRDS(allresu_l, paste0(DEdir,"listDynDE.rds"))
# until here we have padj < 0.1 (not only significant)

# ==========================================================================
# pertinent filter LogFC >= 1  to ** save in table ** for easy interpretation
# yet keeping padj < 0.1
# ==========================================================================
filteredList <- list()
for (n in names(allresu_l)){
  filteredList[[n]] <- allresu_l[[n]] %>% 
    filter(abs(log2FoldChange) >= 1) %>%
    arrange(padj) %>% top_n(500)
}
filteredDF <- dplyr::bind_rows(filteredList)
rownames(genes_df) = genes_df$Geneid
vecsymbols <- sapply(filteredDF$id, function(x) return(genes_df[x,]$symbol) )
filteredDF$symbol = vecsymbols
write.table(filteredDF, file=paste0(DEdir,"filteredDynamicLRT.csv"), 
            sep='\t', col.names=T,row.names=F)
# ==========================================================================
# ====== for spaghetti plot ================

filteredDF <- read.table(paste0(DEdir,"filteredDynamicLRT.csv"),
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
saveRDS(prep,file=paste0(DEdir,"dyn_prep_spaghetti.rds") ) 


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
## end

