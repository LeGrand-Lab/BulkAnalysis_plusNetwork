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
library(cluster)
library(reshape2)
library(RColorBrewer)
library(ggthemes)
library(cowplot)


setwd("~/BulkAnalysis_plusNetwork/")
resdir = "dynamicsIntra_Spec/"

CUTT = 0.5 # the Tau threshold chosen as specific, see 'calc_Tau_Specificity.R'
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
# filter based on Tau index specific or intermediate (>= CUTT):
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

# ========== Venn diagram to show how many shared NOhousekeeping genes:=====
###
venntau <- venn.diagram(x=list("OLD"=tauokgenes.old,
                               "YOUNG"=tauokgenes.young),
                    height=15, width=10, units="cm",
                    col="transparent", fontfamily = "sans",
                    fill=c("violetred", "cornflowerblue"),
                    alpha = 0.5, cex=.8, 
                    cat.col = c("violetred", "blue"),
                    cat.cex = .8, margin=0.1,
                    print.mode=c("raw","percent"), filename=NULL)
          
pdf(paste0(resdir,"tau_shared.pdf"),width=5, height=4)
grid.arrange(gTree(children=venntau), 
       top=paste0("Specific or intermediate genes \n (by Tau >=", CUTT,")"),
       padding=unit(0.2, "line"))
dev.off()
# =========================== end Venn ===================================

# NOTE: GET ONLY ONE 'SUBJECT' BETWEEN OLD OR YOUNG FOR SEPARATE ANALYSIS
# ====================== DE as DYN proxy ===================================
# #
# Function Dynamics by DE (proxy). By default it keeps genes with Tau > 0.3
# by calling funtion 'getNOhousekeeping'
##
customSpecifOnlyDE <- function(here.age, fmat, metadata){
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

#   Dynamics by DE (proxy), use function defined above, then recalculate FDR
# =========================== use custom function ==============================
## ** for OLD **
resu.old <- customSpecifOnlyDE("Old", fmat, metadata)
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

ft.old <- totalOld %>% filter(abs(log2FoldChange)>=2 &
                            repadj <= 0.0005 )
OldIntraDyn_list <- list()
for (i in unique(ft.old$celltype)){
  OldIntraDyn_list[[i]] <- ft.old %>% filter(celltype==i)
}
saveRDS(OldIntraDyn_list, paste0(resdir,"OldIntraDyn_list.rds"))
write.table(ft.old, paste0(resdir,"OldIntraDyn_tab.csv"),
            sep='\t', row.names = T, col.names = T)

## ** for YOUNG **
resu.young <- customSpecifOnlyDE("Young", fmat, metadata)
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
# ============================== end use =============================


# ============================== kmeans  =====================================
# =============== Tuning for optimal 'k' depending on cell type tissue  ==============
here.age = "Old"
needTuning = FALSE # change here if tuning needed
if (needTuning){
  cty <- "ECs" # define here cell type
  testedks = c(2,3,4,5,6,7,8) #testedks = c(4,5,6,7) # for FAPs and sCs use bigger
  agemat = fmat[,str_detect(colnames(fmat), here.age)]
  keep <- apply(agemat, 1, function(row) ifelse(sum(row >=5)>= 3, TRUE, FALSE) )
  agemat <- agemat[keep,]
  here.meta = metadata %>% filter(age==here.age)
  dso <- DESeqDataSetFromMatrix(agemat, here.meta, design=~TYPEtime)
  DYNSPE = readRDS(paste0(resdir, here.age,"IntraDyn_list.rds"))
  subset <- dso[ unique(DYNSPE[[cty]]$id), dso$type==cty ]
  subvsd <- varianceStabilizingTransformation(subset)
  zscores = t(scale(t(subvsd@assays@data@listData[[1]]))) 
  D <- cluster::daisy(zscores) # dissimilarity matrix, to test silhouette
  pdf(paste0(resdir,"silhouette_",cty,".pdf"), width=6, height=6)
  maxsil = c()
  for (Kn in testedks){
    kclusters <- kmeans(zscores, centers=Kn)
    ssi = silhouette(kclusters$cluster, D)
    sussi = summary(ssi)
    plot(silhouette(kclusters$cluster, D), 
         col=1:Kn, border=NA)
    maxsil = c(maxsil,unname(sussi$si.summary[6]) )# the local max
  }
  ggplot(tibble("clusters"=testedks,
                "maxSilhouetteScore"=maxsil), 
         aes(clusters, maxSilhouetteScore)) +
    geom_path() + 
    labs(title=paste(cty,": MaxSilhouette for tested clusters"))
  dev.off()
} else{ print("tuning already done, going to 'clustering' directly")}

# ============================ end tuning   ===============================


# ============================ clustering   ===============================
# from silhouette : 
opticlust <- list("M1"=2,  "M2"=4 ,  "FAPs"=5,  "sCs"=6 , "ECs" =6)
abc = ggthemes::colorblind_pal()(8) # prepare colors for cell types
cellcolors = list(
  "ECs"="#0072B2",
  "FAPs"="#F0E442",
  "M1" = "#D55E00",
  "M2" =  "#CC79A7",
  "Neutro" =  "#009E73",
  "sCs" = "#56B4E9" 
)
here.age = "Old" # "Young
p <- list()
print(paste("dynamicAndSpecific genes found in **", here.age, "** mice"))
DYNSPE = readRDS(paste0(resdir, here.age,"IntraDyn_list.rds"))
for(i in names(DYNSPE)){
  print(paste(i, "==> ", length(unique(DYNSPE[[i]]$id))))
}

ctys <- c("ECs", "FAPs","M1",  "M2" , "sCs")
for (cty in ctys){
  agemat = fmat[,str_detect(colnames(fmat), here.age)]
  keep <- apply(agemat, 1, function(row) ifelse(sum(row >=5)>= 3, TRUE, FALSE) )
  agemat <- agemat[keep,]
  here.meta = metadata %>% filter(age==here.age)
  dso <- DESeqDataSetFromMatrix(agemat, here.meta, design=~TYPEtime)
  subset <- dso[ unique(DYNSPE[[cty]]$id), dso$type==cty ]
  subvsd <- varianceStabilizingTransformation(subset)
  zscores = t(scale(t(subvsd@assays@data@listData[[1]]))) 
  
  kclushere = kmeans(zscores, centers=opticlust[[cty]])
  df_clus = tibble( "id"=names(kclushere$cluster), 
                    "cluster"=paste0("cluster_",kclushere$cluster) )
  write.table(df_clus,paste0(resdir,"clustersKmeans_",here.age,"_",cty,".csv"),
              sep='\t', col.names = T,row.names = F)
  initmeta = here.meta %>% filter(type == cty) %>% select(type,time,TYPEtime,newname)
  if(all(rownames(initmeta)==colnames(zscores))){
    print("rownames(initmeta)==colnames(zscores), ok for reshape2::melt")
    mt = cbind(initmeta,(t(zscores)))
    mtmelt = reshape2::melt(mt, id.vars=c("newname", "time","type","TYPEtime"),
                            variable.name="geneid", value.name="Zscore")  
    mtmelt = full_join(mtmelt, df_clus, by=c("geneid"="id"))
    mtmelt$day = as.numeric(sapply(mtmelt$time, function(x)str_replace(x, "D","")))
    # represent replicate mesures like this ==>  "D4.ENSMUSG00000016529"
    mtmelt$replicate = paste0(mtmelt$time,".",mtmelt$geneid) 
  }else{print("error, zcores colnames and rownames subsetted metadata NOT EQUAL")}
  print(mtmelt[mtmelt$replicate=="D4.ENSMUSG00000016529",])
  # newname time type TYPEtime             geneid     value cluster day             replicate
  # 2 Old.ECs.D4_1   D4  ECs   ECs.D4 ENSMUSG00000016529 0.2461907       6   4 D4.ENSMUSG00000016529
  # 4 Old.ECs.D4_3   D4  ECs   ECs.D4 ENSMUSG00000016529 0.6953550       6   4 D4.ENSMUSG00000016529
  # 5 Old.ECs.D4_2   D4  ECs   ECs.D4 ENSMUSG00000016529 0.6054567       6   4 D4.ENSMUSG00000016529
  lbeller = df_clus %>% group_by(cluster) %>% tally() %>% 
    mutate(clus_nbgenes=paste0(cluster," (n = ",n,")")) %>% pull(clus_nbgenes)
  names(lbeller) = sort(unique(df_clus$cluster))
  
  p[[cty]] <- ggplot(mtmelt, aes(x=day, y=Zscore, group= geneid))  +
    geom_point(size=.1, color="darkgray") + 
    geom_smooth(se=T,
                method="loess", 
                size=0.5,
                alpha=.1, color=cellcolors[[cty]]) + 
    facet_wrap(~cluster, labeller = labeller(cluster=lbeller)) +
    labs(title=paste0("Intra",here.age,", ",cty))
  
}

saveRDS(p, file=paste0(resdir,"complexGGplot_trends",here.age,".rds"))
pdf(paste0(resdir,"SpecificDyn_",here.age,".pdf"), width=7, height=14)
plot_grid(plotlist = p, ncol= 1, rel_heights = c(2,2,1,2,2))
dev.off()


# ================================end code===========================
## Important to know:
# https://support.bioconductor.org/p/114179/ : 
#   In DESeq2, we don't have functionality to test differences between 
#   all consecutive time points, but you can take a look at the stageR
#   package for combining screening p-values with confirmation p-values. 
#   Perhaps this will give you the kind of post-hoc test correction you are looking for.
