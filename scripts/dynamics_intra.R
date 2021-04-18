# compare expression by pairwise between adjacent days
# INSIDE each age separately, as a proxy to "INTRA" dynamics
# Joha GL
##
#library(reticulate)
#use_condaenv(condaenv="~/prog_bio/anaconda3/envs/MYVENV",
#             conda="~/prog_bio/anaconda3/bin/conda")
library(dplyr)
library(tidyverse)
library(VennDiagram)
require(gridExtra) # calling grid.arrange
library("BiocParallel")
register(MulticoreParam(4)) # TODO:  set n of cores depending of available
library(DESeq2)

setwd("~/BulkAnalysis_plusNetwork/")
resdir = "dynamicsIntra/"

prefil_cou <- "data/prefiltered_counts.rds"
metadata.rds <- "data/metadata.rds"
genes_df <- read.table("data/genesinfo.csv",sep="\t",header=T)

fmat <- readRDS(prefil_cou)
metadata <- readRDS(metadata.rds)
metadata <- metadata %>% mutate(TYPEtime=paste0(type,".",time)) 

# rows to keep
keep <- apply(fmat, 1, function(row) ifelse(sum(row >=5)>= 3, TRUE, FALSE) )
fmat <- fmat[keep,]

# NOTE: we are interested only in most "specific" genes dynamics, 
# filter based on Tau index specific or intermediate (> 0.3):
# THIS WILL DISCARD HOUSEKEEPING GENES
getNOhousekeeping <- function(age){
  tauokgenes <- c()
  myregex = paste0(age,".*\\.txt$")
  for (i in list.files("Tau/", pattern=myregex)){
    itext <- read.table(paste0("Tau/",i), sep='\t', header=T, row.names = 1)
    tmpvec <- itext %>% filter(class !='housekeeping') %>% pull(id)
    tauokgenes <- c(tauokgenes, tmpvec)
  }
  return(unique(tauokgenes))
}

tauokgenes.old <- getNOhousekeeping("Old")
tauokgenes.young <- getNOhousekeeping("Young")

# Venn diagram to show how many shared NOhousekeeping genes:
venntau <- venn.diagram(x=list("OLD"=tauokgenes.old,
                    "YOUNG"=tauokgenes.young),
                    height=15, width=10, units="cm",
                    col="transparent", fontfamily = "sans",
                    fill=c("violetred", "cornflowerblue"),
                    alpha = 0.5, cex=.8, 
                    cat.col = c("violetred", "blue"),
                    cat.cex = .8, margin=0.1,
                    print.mode=c("raw","percent"), filename=NULL)
          
pdf(paste0(resdir,"tau_shared_correct.pdf"),width=5, height=4)
grid.arrange(gTree(children=venntau), 
       top="Specific or intermediate genes \n (by Tau > 0.3)",
       padding=unit(0.2, "line"))
dev.off()

# NOTE: GET ONLY ONE 'SUBJECT' BETWEEN OLD OR YOUNG FOR SEPARATE ANALYSIS
# ===================================================================
# #
# Function Dynamics by DE (proxy). By default it keeps genes with Tau > 0.3
# by calling funtion 'getNOhousekeeping'
##
customLongDE <- function(here.age, fmat, metadata){
  agemat = fmat[,str_detect(colnames(fmat), here.age)]
  keep <- apply(agemat, 1, function(row) ifelse(sum(row >=5)>= 3, TRUE, FALSE) )
  agemat <- agemat[keep,]
  # filter matrix depending on  tissue genes specificity (Tau ! ) :
  here.tau.genes = getNOhousekeeping(here.age)
  agemat = agemat[rownames(agemat) %in% here.tau.genes,]
  here.meta = metadata %>% filter(age==here.age)
  dso <- DESeqDataSetFromMatrix(agemat, here.meta, design=~TYPEtime)
  d <- DESeq2::estimateSizeFactors(dso,type="ratio")
  d <- DESeq2::estimateDispersions(d)
  resulist = list()
  cts = unique(metadata$type)
  # to make contrasts proper to each celltype, initialize empty  matrix design
  MyMat <- array(0, dim=c(dim(here.meta)[1], length(unique(here.meta$TYPEtime))))
  rownames(MyMat) <- rownames(here.meta)
  colnames(MyMat) <- unique(here.meta$TYPEtime)                
  for (icts in cts){
    print(paste("    ===> current celltype :", icts ))
    tpoints = here.meta %>% filter(type == icts) %>% pull(time) %>% unique()
    tpoints = sort(tpoints)
    if (length(tpoints) > 1){
      for (i in 2:length(tpoints)){
        iref = paste0(icts, ".", tpoints[i-1])  # an example is 'M2.D0'
        inoref = paste0(icts, ".", tpoints[i])
        # *!* dynamically get matrix design *!*
        i.mat.mod <- MyMat  
        i.mat.mod[str_detect(rownames(i.mat.mod), iref) , iref  ] <- 1
        i.mat.mod[str_detect(rownames(i.mat.mod), inoref) , inoref  ] <- 1
        i.mat.mod <- i.mat.mod[,colSums(abs(i.mat.mod)) > 0]
        #intercept <- apply(i.mat.mod, 1, sum)
        #imat.mod <- as.array(cbind(intercept,i.mat.mod))     
        d1 <- nbinomWaldTest(d, modelMatrix = i.mat.mod) 
        rd1 <- results(d1, contrast = list(inoref,iref), parallel=T)
        rd1$id <- rownames(rd1)
        rd1$contrast <- paste0(c(inoref,iref),collapse="_vs_")
        rd1$celltype <- icts
        resulist[[paste0(c(inoref,iref),collapse="_vs_")]] <- as_tibble(rd1)  
        # NOTE: not yet converted ensembl id to symbol
      }# end for
    }else{ print("contrasts NOT feasible, only has a unique timePoint")}
  }
  return(resulist)
}  # end function

#   Dynamics by DE (proxy)
# ==========================================================================
## for OLD
resu.old <- customLongDE("Old", fmat, metadata)
totalOld <- dplyr::bind_rows(resu.old)
dim(totalOld) # 95988 9  # this means the total of tests, so, MUST:
# recalculate FDR, because multiple testing performed above)
pvals = totalOld$pvalue
totalOld$repadj = p.adjust(pvals, method = "BH", n = length(pvals))
# plot this FDR improvement : 
pdf(paste0(resdir,"recalcFDRintra.pdf"), height=6, width= 6)
ggplot(totalOld, aes(x= pvalue)) +
  geom_histogram(fill="gold", alpha=0.5) +
  geom_histogram(data=totalOld, aes(x=padj), fill="red", alpha=0.4) +
  geom_histogram(data=totalOld, aes(x=repadj), fill="green", alpha=0.3) +
  facet_wrap(~contrast) + 
  labs(title = "Recalculated FDR across contrasts",
        x = "pvalue(gold), padj(red), re-padj(green)")
dev.off()
# filterdataframe by lfc >= 2, repadj < 0.005
ft.old <- totalOld %>% filter(abs(log2FoldChange)>=2 &
                            repadj <= 0.0005 )
OldIntraDyn_list <- list()
for (i in unique(ft.old$celltype)){
  OldIntraDyn_list[[i]] <- ft.old %>% filter(celltype==i)
}
saveRDS(OldIntraDyn_list, paste0(resdir,"OldIntraDyn_list.rds"))
write.table(ft.old, paste0(resdir,"OldIntraDyn_tab.csv"),
            sep='\t', row.names = T, col.names = T)

## for YOUNG
resu.young <- customLongDE("Young", fmat, metadata)
totalYoung <- dplyr::bind_rows(resu.young)
dim(totalYoung) 
# recalculate FDR, because multiple testing performed above
pvals.young = totalYoung$pvalue
totalYoung$repadj = p.adjust(pvals.young, method = "BH", n = length(pvals.young))
# filterdataframe by lfc >= 2, repadj < 0.005
ft.young <- totalYoung %>% filter(abs(log2FoldChange)>=2 & repadj <= 0.0005 )
YoungIntraDyn_list <- list()
for (i in unique(ft.young$celltype)){
  YoungIntraDyn_list[[i]] <- ft.young %>% filter(celltype==i)
}
saveRDS(YoungIntraDyn_list, paste0(resdir,"YoungIntraDyn_list.rds"))
write.table(ft.young, paste0(resdir,"YoungIntraDyn_tab.csv"),
            sep='\t', row.names = T, col.names = T)
# ==========================================================================
# end finding most dyn genes by proxy

# ======================== kmeans CLUSTERING ===============================
# ==========================================================================
# first testing OLD SUBJECT
age = "Old"
vstmat <- get vst stuff !!!entire matrix # old
# calculate z score matrix , as ideal for clustering:
zmat = t(scale(t(vstmat)))
# discard all that is not dynamic
zmat <- zmat[onlydyn,]
# split into celltype sub matrices list
ctlist <- list()
clustres <- list()
for (ict in cts){
  ctlist[[ict]] <- zmat[,str_detect(colnames(zmat),ict)]
  clustres <- 
}


## Important to know:
# https://support.bioconductor.org/p/114179/ : 
#   In DESeq2, we don't have functionality to test differences between 
#   all consecutive time points, but you can take a look at the stageR
#   package for combining screening p-values with confirmation p-values. 
#   Perhaps this will give you the kind of post-hoc test correction you are looking for.


