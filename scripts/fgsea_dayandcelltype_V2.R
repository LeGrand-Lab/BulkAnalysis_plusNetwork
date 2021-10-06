#  Runs GSEA by day and celltype : VERSION 2 
# from 'static' DEGs
# Saves rds objects into exam_INTER_conditions/static/GSEA/rds/
# * Major modification : only REACTOME pathways tested *
# note: filtered rds file contains only top 15 up and top 15 down pathways
# --
# johaGL
library(fgsea)
library(tidyverse)
library(msigdbr)
library(cowplot)
setwd("~/BulkAnalysis_plusNetwork/")

daysv = c('D0', 'D2', 'D4', 'D7')
DEfile = "rds/shot_rds_full.rds"
odir = "exam_INTER_conditions/static/"
infoinputfile = "GSEA/csv/infoinput-V2"  # csv or md_txt extensions 
pdfplotfile = "preGSEA_byday_bycelltype-V2.pdf"
gseards = "GSEA/rds/"
gseaoutfull = "fgsea_bd_bct_full-V2.rds" 
gseaoutfiltered = "GSEA/rds/fgsea_bd_bct_filtered-V2.rds"
consensusfile <- "Tau/conseTau_ensemblid.rds"


consensustau <- readRDS(consensusfile)
DEdf <- readRDS(paste0(odir, DEfile))
summary(DEdf)
# names(DEdf)
# [1] "baseMean"       "log2FoldChange" "lfcSE"          "pvalue"         "padj"           "id"             "day"           
# [8] "type"           "symbol"        
# > names(consensustau)
# [1] "D0" "D2" "D4" "D7"

if (!"symbol" %in% colnames(DEdf)){
  genes_df <- read.table("data/genesinfo.csv",sep="\t",header=T)
  DEdf$symbol <- genes_df[match(DEdf$id, genes_df$Geneid),]$symbol
  head(DEdf)
  print("added symbols to dataframe")
}else{ print("symbol already in dataframe")}

concatenateDEwithTau <- function(DEdf, consensustau){
  mixedfull <- list()
  for (d in names(consensustau)){
    print(paste(" %%%%%%%%%%%%%", d, "%%%%%%%%%%%%%"))
    dem.here <- DEdf %>% filter(day == d) %>% filter(!is.na(padj))
    tau.here <- consensustau[[d]]
    tmp_ <- list()
    for (c in unique(dem.here$type)){
      print(paste("     ", c))
      mix <- full_join(
          dem.here %>% filter(type==c),
          tau.here %>% filter(whichMAX==c),
          by = "id"  )# end _join
      tmp_[[c]] <- mix
    } # end for
    mixedfull[[d]] <- bind_rows(tmp_)
  } # end for
  return(mixedfull)
}
  
mixedfull <- concatenateDEwithTau(DEdf, consensustau)
# tail(mixedfull$D2 %>% filter(!is.na(Tau) & !is.na(padj)) %>% arrange(desc(padj)))
# baseMean log2FoldChange lfcSE   pvalue     padj id                 day   type   symbol.x symbol.y Tau               class        whichMAX nbMAX exclusiveOld
# <dbl>          <dbl> <dbl>    <dbl>    <dbl> <chr>              <chr> <chr>  <chr>    <chr>    <chr>             <chr>        <chr>    <chr> <chr>       
#   1     429.           2.02 0.221 6.49e-21 1.35e-17 ENSMUSG00000026768 D2    sCs    Itga8    Itga8    0.767054165861923 intermediate sCs      1     0  





# initial plots

givemeggplot <- function(mix){
  tx <- mix %>% mutate(Tauplo = ifelse(is.na(Tau), -1, Tau)) %>% filter(!is.na(padj))
tx$Tauplo <- as.numeric(tx$Tauplo)
ggplot(tx, aes(x=-padj, y=log2FoldChange, color=Tauplo )) + 
  geom_point(size=.5) 
ggplot(tx, aes(x=Tauplo, y=log2FoldChange, color=padj<=0.1 )) + 
  geom_point(size=.5) + scale_color_brewer(palette='Dark2') + facet_grid(~padj <0.05)}

givemeggplot


# Mosaicoooooooooooooooooooooooooooooooooooooooooo




# === split by day and keep in a list of dataframes (each df ~1000 genes): =====
DE_l = list()  #  
for (k in daysv){
  tmp <- DEdf %>% filter(day == k)  %>% 
    mutate(sens= ifelse(log2FoldChange < 0,"down", "up"))
  # keep top genes by abslfc and pval
  tmp_up <- tmp %>% filter(sens == "up") %>% group_by(type) %>%
    arrange(padj, .by_group = TRUE) %>% slice_min(padj, n=M/2, with_ties = F)
  tmp_down <- tmp %>% filter(sens == "down") %>% group_by(type)   %>%
    arrange(padj,  .by_group = TRUE) %>% slice_min(padj, n=M/2, with_ties = F)
  DE_l[[k]] <-  rbind(tmp_up, tmp_down)
}
infoinput <- data.frame("day_celltype" = c(), "inputsize"=as.integer(c()), 
                        "maxpadj" = c(), "minabslfc" = c())
for (k in daysv){
  print("")
  print(paste(":::::",k,":::::"))
  cts <- unique(DE_l[[k]]$type) 
  tmpplots_ <- list()
  for (CT in cts) {
    print(paste("   --->", CT))
    here.df <- DE_l[[k]] %>% filter(type == CT) %>% arrange(desc(log2FoldChange)) %>%
      select(log2FoldChange, symbol, padj)
    nbg = dim(here.df)[1]
    maxp = max(here.df$padj)
    minalf = min(abs(here.df$log2FoldChange))
    infoinput <- rbind(infoinput, c(paste0(k,"_",CT),
                                    nbg, maxp, minalf))
  } }
colnames(infoinput) <- c("day_celltype" , "inputsize", 
                         "maxpadj" , "minabslfc" )
infoinput$inputsize <- as.integer(infoinput$inputsize)
write.table(infoinput, paste0(odir, infoinputfile, ".csv"), sep='\t', 
            col.names = T, row.names = F)
write.table(infoinput %>% mutate(day_celltype = paste("|  ", day_celltype)) %>%
              mutate(minabslfc = paste(minabslfc, " |")), 
            file = paste0(odir, infoinputfile, "_md.txt"), sep=" | ", 
            col.names = T, row.names = F)

# ======================== perform Gsea on top ranked genes ====================

# visualize and save barplots before running fgsea
plots_ <- list()
for (k in daysv){
  cts <- unique(DE_l[[k]]$type) 
  tmpplots_ <- list()
  for (CT in cts) {
    here.df <- DE_l[[k]] %>% filter(type == CT) %>% arrange(desc(log2FoldChange)) %>%
      select(log2FoldChange, symbol)
    gseagenes = here.df %>% pull(log2FoldChange)
    names(gseagenes) = here.df$symbol
    barplot(gseagenes, main=paste(k, CT))
    # same barplot but on ggplot: 
    here.df <- here.df %>% mutate(seq=1:n())
    new <- c(here.df[1,]$symbol) # new labels, first element
    saut = 40
    addsum = 40
    for (i in 2:dim(here.df)[1]){
      if( i == saut){
        new = c(new, here.df[i,]$symbol )
        saut = saut + addsum
      }else{ new = c(new, "") } } # end for
    gplo <- ggplot(here.df , aes(x=reorder(symbol, log2FoldChange),y=log2FoldChange)) +
      geom_col() + 
      theme( axis.ticks.x = element_blank(),
             axis.text.x = element_text(angle = 90)) +
      scale_x_discrete(labels= new) + 
      labs(title=paste(k, CT), x = "genes") 
    tmpplots_[[CT]] <- gplo
  }
  plots_[[k]] <- tmpplots_
}
pdf(paste0(odir, "GSEA/", pdfplotfile), width=12, height=30)
plots_[['D0']][["e1"]] <- NULL
plot_grid(
  plot_grid(plotlist = plots_[['D0']], nrow=2, ncol=3),
  plot_grid(plotlist = plots_[['D2']], nrow= 2, ncol = 3),
  plot_grid(NULL),
  plot_grid(plotlist = plots_[['D4']], nrow=2, ncol = 3),
  plot_grid(NULL),
  plot_grid(plotlist = plots_[['D7']], nrow=2, ncol = 3),
  nrow = 6, rel_heights = c(3,3,1,3,1,3)
)
dev.off()

# ======================== RUNNING GSEA =======================================
thegmt <- msigdbr(species = "Mus musculus", 
                  category = 'C2', 
                  subcategory=c('CP:REACTOME'))
msigdbr_list = split(x = thegmt$gene_symbol, f = thegmt$gs_name)

pathsFull_l <- list("D0"=list(),"D2"=list(),"D4"=list(),"D7"=list())
pathsFiltered  <- list("D0"=list(),"D2"=list(),"D4"=list(),"D7"=list())
for (k in daysv){
  cts <- unique(DE_l[[k]]$type) 
  dd <- list()
  dd_f <- list()
  for (CT in cts){
    set.seed(42)
    here.df <- DE_l[[k]] %>% filter(type == CT) %>% arrange(desc(log2FoldChange)) %>%
      select(log2FoldChange, symbol)
    print(c(k, CT))
    gseagenes = here.df %>% pull(log2FoldChange)
    names(gseagenes) <- here.df$symbol
    pathsFull_l[[k]] = list()
    pathsFiltered[[k]] = list()
    print("running fgseaMultilevel, nperm not needed")
    fgseaRes <- fgsea::fgsea(pathways = msigdbr_list, 
                             stats = gseagenes,
                             minSize=3,
                             maxSize=Inf) 
    topPathwaysUp <- fgseaRes[ES>0][head(order(padj),n=15), ]
    topPathwaysDown <- fgseaRes[ES<0][head(order(padj), n=15), ]
    combipath <- rbind(topPathwaysUp, topPathwaysDown)
    combipath <- combipath %>% 
      mutate(sens = ifelse(ES > 0, "UP", "DOWN" ))
    dd[[CT]] <- fgseaRes
    dd_f[[CT]] <- combipath
  }
  pathsFull_l[[k]] <- dd
  pathsFiltered[[k]] <- dd_f
}

saveRDS(pathsFull_l, file=paste0(odir,gseards, gseaoutfull ))
saveRDS(pathsFiltered, file=paste0(odir,gseards, gseaoutfiltered  ))


# ================== printing GSEA min padj results ==========================
XX = readRDS(paste0(odir,gseardsl, gseaoutfull))
padjofpaths = array(NA, dim=c(4,6))
rownames(padjofpaths) = daysv
colnames(padjofpaths) = names(XX[['D2']])
for (d in daysv){
  for (n in names(XX[[d]])){
    print(n)
    print(min(XX[[d]][[n]]$padj))
    padjofpaths[d, n] <- min(XX[[d]][[n]]$padj)
  }
}
# =========================================================================
