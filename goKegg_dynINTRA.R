# joha GL 2021
library(dplyr)
library(tidyverse)
library(ggplot2)
library(cowplot)
library(forcats)
library(RColorBrewer)
library(MASS)
library(pheatmap)
library(ggsci) # publishing palettes
library(cowplot)
library(gridExtra)
library(reshape2)
library(ggthemes)
library(scales) # 'viridis_pal' et autres pal
library(ggrepel) # for labels to points
library(org.Mm.eg.db)  
library(AnnotationDbi)
library(edgeR)


setwd("~/BulkAnalysis_plusNetwork")
resdir = "plotsDE/"
DEdir = "signaturestypes/"
prefil_cou <- "data/prefiltered_counts.rds"
metadata.rds <- "data/metadata.rds"
# https://cran.r-project.org/web/packages/msigdbr/vignettes/msigdbr-intro.html

vrs = "vB"  # TODO: file version 

ens2entrez <- function(vectorensemblids){
  res <- AnnotationDbi::select(org.Mm.eg.db,
                               key=vectorensemblids,
                               columns="ENTREZID", 
                               keytype = "ENSEMBL")
  return(res) # a dataframe columns ENSEMBL and ENTREZID
}
entrez2symbol <- function(vectorensemblids){
  res <- AnnotationDbi::select(org.Mm.eg.db,
                               key=vectorensemblids,
                               columns="SYMBOL", 
                               keytype = "ENTREZID")
  return(res) # a dataframe columns ENSEMBL and ENTREZID
}
# ========================  part Edger with limma kegg and go

fmat <- readRDS(prefil_cou)
metadata <- readRDS(metadata.rds)
metadata <- metadata %>% mutate(timetype=paste0(time,".",type)) 
# rows to keep
keep <- apply(fmat, 1, function(row) ifelse(sum(row >=5)>= 3, TRUE, FALSE) )
fmat <- fmat[keep,]
allct <- sort(unique(metadata$type))
allct <- allct[allct !="Neutro"]  #  "ECs"  "FAPs" "M1"   "M2"   "sCs" 
ages=c("Young","Old")

# kegg by pairwise: 
for (ag in ages){
  resu_l <- list()
  gokegg <- list()
  for (ct in allct){
    agectmx <- fmat[,str_detect(colnames(fmat), ag) & str_detect(colnames(fmat), ct)]
    acmeta <- metadata %>% filter(age==ag & type==ct)
    timeps <- sort(unique(acmeta$time))
    Ntp = length(timeps)   #  this variable is used downstream in this loop
    print(timeps)
    timepoints <- factor(acmeta[match(colnames(agectmx),acmeta$newname),]$time,
                         levels=timeps)
    # avoid multiple ensembl to 1 common entrez :
    ens_entrez <- ens2entrez(rownames(agectmx)) %>% as_tibble() %>%
      group_by(ENTREZID) %>% filter(!(n_distinct(ENSEMBL) >= 2))
    # delete rows not being maped to entrezid
    agectmx <- agectmx[ens_entrez$ENSEMBL,]
    rownames(agectmx) <- ens_entrez[match(rownames(agectmx),
                                          ens_entrez$ENSEMBL),]$ENTREZID
    y <- DGEList(counts=agectmx, group=timepoints)
    keep <- filterByExpr(y)
    y <- y[keep, , keep.lib.sizes=FALSE]
    y <- calcNormFactors(y)
    basecontrast <- rep(0,Ntp)
    design <- model.matrix(~0+timepoints)
    y <- estimateDisp(y,design)
    fit <- glmQLFit(y, design)
    tmp_l <- list()
    go_l <- list()
    ke_l <- list()
    for (k in 2:Ntp){
      j <- k-1
      basecontrast[j] <- -1
      basecontrast[k] <- 1
      print(basecontrast)
      qlf <- glmQLFTest(fit, contrast=basecontrast)
      dynamicN <- ifelse(Ntp==2, 300 ,ifelse(Ntp==3, 150, 80))
      keg <- kegga(qlf,species="Mm")
      go <- goana(qlf, species="Mm")
      gr <- topGO(go, ontology=c("BP","CC", "MF"), 
                   number=dynamicN)
      kr <- topKEGG(keg, number=Inf)
      tr = as.data.frame(topTags(qlf, sort.by="logFC",p.value=0.3, n=dynamicN))
      tr$type = ct
      thiscontrast = paste0(c(timeps[k],timeps[j]),collapse="vs")
      tr$contrast = thiscontrast
      tr$id = rownames(tr)
      tmp_l[[thiscontrast]] <- tibble(tr)
      go_l[[thiscontrast]] <- gr
      ke_l[[thiscontrast]] <- kr
      basecontrast <- rep(0,Ntp) # reinitialize base contrast
    }
    gokegg[[paste0(ct,"_go")]] <- go_l
    gokegg[[paste0(ct,"_kegg")]] <- ke_l
    resu_l[[ct]] <- bind_rows(tmp_l)
  }#end for ct
  kkk <- bind_rows(resu_l)
  write.table(kkk, file=paste0(DEdir,"edger_dynINTRA_",ag,vrs,".txt"), sep='\t',
              col.names=T, row.names=F)
  saveRDS(gokegg, file=paste0(DEdir,"edger_dynINTRAkeggo",ag,vrs,".rds"))
}#end for


# ========================  doing in form of heat map, only for M1

allkeggenes = getGeneKEGGLinks(species.KEGG = "mmu", convert = FALSE)
dim(allkeggenes) # 36015     2
gokegg <- readRDS( file=paste0(DEdir,"edger_dynINTRAkeggo","Old",vrs,".rds"))
DYNgenes <- read.table(paste0(DEdir,"edger_dynINTRA_","Old",vrs,".txt"), sep='\t',
                       header=T)
DYNkegg_df <- allkeggenes[allkeggenes$GeneID %in% unique(DYNgenes$id),]
dim(DYNkegg_df) # [1] 2032 2  (only including our more dynamic genes )

names(gokegg[["M1_kegg"]])
trast <- "D4vsD2"
tmpkegg_t_t <- gokegg[["M1_kegg"]][[trast]] %>% 
  filter(N >= mean(tmpkegg_t_t$N))  #  has unique pathway ids as rownames ('path:mmu04933')
localgenes <- DYNgenes %>% filter(type=="M1" & contrast==trast)
DYNkegg_t <- allkeggenes %>% filter(PathwayID %in% rownames(tmpkegg_t_t) &
            GeneID %in% unique(localgenes$id))

matrixlfc <- array(0,dim=c(length(unique(localgenes$id)), 
                           length(unique(DYNkegg_t$PathwayID))  ) )
rownames(matrixlfc) <- unique(localgenes$id)
colnames(matrixlfc) <- unique(DYNkegg_t$PathwayID)
for (gene in localgenes$id){
  #print(localgenes[localgenes$id==gene,]$logFC)
  thesepaths <- DYNkegg_t[DYNkegg_t$GeneID == gene,] 
  for (path in thesepaths$PathwayID){
    matrixlfc[as.character(gene),path] <- localgenes[localgenes$id==gene,]$logFC
  }
}
# transform column and rows to names useful to eyes:
labeledm <- matrixlfc
colnames(labeledm) <- tmpkegg_t_t[colnames(matrixlfc),]$Pathway
rownames(labeledm) <- entrez2symbol(rownames(matrixlfc))$SYMBOL
keep <- apply(labeledm,1, function(x) sum(x)==0)
names(keep) <- as.character(names(keep))

colorgenes <- rep(c("red","blue"), length(rownames(labeledm))+1)

pheatmap(labeledm[!keep,], border_color="white",
         clustering_distance_cols = "correlation",
         color = colorRampPalette(rev(brewer.pal(n = 7, 
                                                 name = "PuOr")))(100),
         fontsize = 7,
         angle_col = 45,
         main = "M1 : D4vsD2")


p <- pheatmap(labeledm[!keep,], border_color="white",
              clustering_distance_cols = "correlation",
              color = colorRampPalette(rev(brewer.pal(n = 7, 
                                                      name = "PuOr")))(100),
              fontsize = 7,
              angle_col = 45,
              main = "M1 : D4vsD2") + theme_bw()

my_gtable = p$gtable
my_gtable$grobs[[6]]$gp=grid::gpar(col=colorgenes, fontsize=8)
pdf(paste0(resdir,"dynamicTOPgenes_M1.pdf"))
plot(my_gtable)
dev.off()

# ========================  doing in form of Alluvial

library(ggalluvial) ## include in pkg install 
library(purrr)  ##include in pkg install 

options(stringsAsFactors = F)

#Generate dummy data
myData = map_df(paste0("Person", 1:9), function(x){
  data.frame(     Person = x, 
             Characteristic = paste0("Characteristic ",
                                     sample(1:5, sample(2:5, 1))))
})
myData$other <- rep(c("otherC","otherB","otherA"), (dim(myData)[1]/3) )
myData$lfc <- rep(c(-2,-1,0,2,3,1.5, 2.2), ((dim(myData)[1]/6)))[1:33]
 #https://cran.r-project.org/web/packages/ggalluvial/vignettes/labels.html
#Create the alluvial plot
ggplot(data = myData,
       aes(axis1 = other, axis2=Person, axis3 = Characteristic, y = 0.5)) +
  scale_x_discrete(expand = c(.4, 0)) +
  geom_alluvium(aes(fill = lfc), show.legend = FALSE) +
  geom_stratum() + geom_text(aes(label = Person),
                             stat = "stratum", size = 3) +
  theme_void() + ggtitle("Alluvial plot") + 
  scale_fill_distiller(type = 'seq', palette = 'GnBu',
                       direction = -1, 
                       labels = comma, 
                       guide = "colourbar", 
                       aesthetics = "fill") + 
  ggrepel::geom_text_repel(
    aes(label =  as.character(other)),
    stat = "stratum", size = 4, direction = "y", nudge_x = -.5
  ) + 
  ggrepel::geom_text_repel(
    aes(label =  as.character(Characteristic)),
    stat = "stratum", size = 4, direction = "y", nudge_x = .5
  )
  #guides(fill = guide_colorbar(barwidth = 0.5, barheight = 10, title="Value"))
  
  



ggplot(myData,
       aes(x = other, stratum = Person, alluvium = Characteristic, y = 0.5,
           fill = response)) +
  scale_x_discrete(expand = c(.4, 0)) +
  geom_flow(width = 1/4) +
  geom_stratum(alpha = .5, width = 1/4) +
  scale_linetype_manual(values = c("blank", "solid")) +
  ggrepel::geom_text_repel(
    aes(label = ifelse(as.numeric(survey) == 1, as.character(response), NA)),
    stat = "stratum", size = 4, direction = "y", nudge_x = -.5
  ) +
  ggrepel::geom_text_repel(
    aes(label = ifelse(as.numeric(survey) == 3, as.character(response), NA)),
    stat = "stratum", size = 4, direction = "y", nudge_x = .5
  ) +
  theme(legend.position = "none") +
  ggtitle("vaccination survey responses", "labeled using `geom_text_repel()`")


colorRampPalette(rev(brewer.pal(n = 7, 
                                name = "PuOr")))(100),
#  Note: VERY USEFUL FOR GSEA:
# library("msigdbr") # install.packages("msigdbr")
#all_gene_sets = msigdbr(species = "Mus musculus")
#View(msigdbr_collections())
#setsli <- list()
#setsli[["H"]] = msigdbr(species = "Mus musculus", category = "H")
#setsli[["C2"]] = msigdbr(species = "Mus musculus", category = "C2")
#saveRDS(setsli, file=paste0(DEdir,"DATABASE.rds"))
#setsli <- readRDS(paste0(DEdir,"DATABASE.rds"))
#View(as_tibble(unique(setsli[["C2"]]$gs_name))) # str_detec REACTOME and KEGG to extract
#setsli[["H"]] is hallmark . REACTOME and KEGG could be enough