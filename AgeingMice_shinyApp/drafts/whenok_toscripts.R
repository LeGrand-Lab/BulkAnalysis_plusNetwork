#  TODO: retrieve pvalue
# --
# johaGL
library(fgsea)
setwd("~/BulkAnalysis_plusNetwork/")
# take top 200 ranked genes by pvalue and abs(l2fc))
# and perform gsea like analysis with gprofiler2

odir = "exam_INTER_conditions/static/"
shot_t = "shot_dataframe_unfiltered.csv"   # this is soft filtered version
# shot_t_f = "shot_dataframe.csv"   # this is filtered version

DEuf <- read.table(paste0(odir, shot_t), header=T, sep='\t')
# DEfi <- read.table(paste0(odir, shot_t_f), header=T, sep=',')
summary(DEuf)
# gene symbol missing, add it
genes_df <- read.table("data/genesinfo.csv",sep="\t",header=T)

DEuf$symbol <- genes_df[match(DEuf$id, genes_df$Geneid),]$symbol
head(DEuf)

# split by day and keep in a list that will be our query:
DE_l = list()  #  
for (k in c('D0','D2', 'D4', 'D7')){
  tmp <- DEuf %>% filter(day == k)  %>% 
    mutate(absLFC = abs(log2FoldChange)) %>% mutate(sens= ifelse(log2FoldChange < 0,"down", "up"))
  # keep top genes by abslfc and pval
  tmp_up <- tmp %>% filter(sens == "up") %>% group_by(type)   %>%
    arrange(padj, desc(absLFC), .by_group = TRUE) %>% slice_min(padj, n=50)
  tmp_down <- tmp %>% filter(sens == "down") %>% group_by(type)   %>%
    arrange(padj, desc(absLFC), .by_group = TRUE) %>% slice_min(padj, n=50)
  DE_l[[k]] <-  rbind(tmp_up,tmp_down)
}

# ======================== perform Gsea on top ranked genes
thegmt <- read.table(paste0("stock_gmtfiles/","Hallmark_React.gmt"), sep='\t',
                     header=T)
msigdbr_list = split(x = thegmt$gene_symbol, f = thegmt$gs_name)

set.seed(42)

pathsFull_l <- list("D0"=list(),"D2"=list(),"D4"=list(),"D7"=list())
pathsFiltered  <- list("D0"=list(),"D2"=list(),"D4"=list(),"D7"=list())
for (k in c('D0', 'D2', 'D4', 'D7')){
  cts <- unique(DE_l[[k]]$type) 
  dd <- list()
  dd_f <- list()
  for (CT in cts){
    here.df <- DE_l[[k]] %>% filter(type == CT)
    print(c(k, CT))
    gseagenes = here.df %>% arrange(desc(absLFC)) %>% pull(log2FoldChange)
    names(gseagenes) <- here.df$symbol
    print(gseagenes)
    barplot(sort(gseagenes))
    pathsFull_l[[k]] = list()
    pathsFiltered[[k]] = list()
    fgseaRes <- fgsea::fgsea(pathways = msigdbr_list, 
                             stats = gseagenes,
                             minSize=3,
                             maxSize=Inf, nperm = 100000) 
    topPathwaysUp <- fgseaRes[ES>0][head(order(padj),n=15), ]
    topPathwaysDown <- fgseaRes[ES<0][head(order(padj), n=15), ]
    combipath <- rbind(topPathwaysUp, topPathwaysDown)
    combipath <- combipath %>% 
      mutate(path4graph = str_replace(pathway, "REACTOME_", "")) %>% 
      mutate(sens = ifelse(ES > 0, "UP", "DOWN" ))
    dd[[CT]] <- fgseaRes
    dd_f[[CT]] <- combipath
  }
  pathsFull_l[[k]] <- dd
  pathsFiltered[[k]] <- dd_f
}

saveRDS(pathsFull_l, file="AgeingMice_shinyApp/drafts/fgseaByDay_full.rds" )
saveRDS(pathsFiltered, file="AgeingMice_shinyApp/drafts/fgseaByDay_filtered.rds" )

pathsFiltered = readRDS("AgeingMice_shinyApp/drafts/fgseaByDay_filtered.rds" )

## prepare matrices for heatmaps
mhp = list()
for (k in c('D0','D2', 'D4', 'D7')){
  cts <- unique(DE_l[[k]]$type) 
  mhp[[k]] <- list()
  for (CT in cts){
    mhp[[k]][[CT]] <- list()
    for(S in c("UP","DOWN")){
      here.data = as.data.frame(pathsFiltered[[k]][[CT]])
      here.data <- here.data %>% filter(sens == 'UP')
      here.dico = list()
      if(length(here.data$path4graph) == length(unique(here.data$path4graph) )){
        print("yjjsjjsjj")
        rownames(here.data) =  here.data$path4graph
        aggreg_genes = c()
        for (path in here.data$path4graph){
          tmpstr <- str_replace_all(here.data[path,]$leadingEdge, 
                                    c('c\\("' = '',   '"' = '', '\\)' = ''))
          here.dico[[path]] <- unname(unlist(str_split(tmpstr,", ")))
          aggreg_genes <- c(aggreg_genes, here.dico[[path]] )
        }
        colnamesgenes = unique(aggreg_genes)
        matrx = matrix(NA, length(here.data$path4graph), length(colnamesgenes))
        rownames(matrx) = here.data$path4graph
        colnames(matrx) = colnamesgenes
        for (path in here.data$path4graph){
          localgenes <- here.dico[[path]]
          for (g in localgenes){
            matrx[path,g] <-  here.data[path,]$NES # if more than once, yields 2 or more
          }
        }
        # clear rows and columns containing 2 or less values
        matrx <- matrx[rownames(matrx) != "DISEASE",]          # take away "DISEASE" is stupid
        rtokeep <- apply(matrx, 1, function(x) sum(!is.na(x)) >= 4)
        ctokeep <- apply(matrx, 2, function(x) sum(!is.na(x)) >= 4)
        newm <- matrx[rtokeep,ctokeep]
        
        mhp[[k]][[CT]][[S]] <- newm
      } else {
        print(paste("error, this", k, CT ,"  contains repeated pathwaynames") )
        stop()
      } # endif
      
    } # end for sens UP or DOWN
    
  }}


View(mhp[["D7"]][["M2"]][["UP"]])

saveRDS(pathsFiltered, file="AgeingMice_shinyApp/drafts/matrices4heatmaps.rds" )



# ## initial test fgsea:
# set.seed(42)
# pathsFull_l = list()
# pathsFiltered = list()
# k = 'D0'
# cts <- unique(DE_l[[k]]$type) 
# CT <- 'ECs'
# here.df <- DE_l[[k]] %>% filter(type == CT)
# gseagenes = here.df %>% arrange(desc(absLFC)) %>% pull(log2FoldChange)
# names(gseagenes) <- here.df$symbol
# barplot(sort(gseagenes))
# pathsFull_l[[k]] = list()
# pathsFiltered[[k]] = list()
# fgseaRes <- fgsea::fgsea(pathways = msigdbr_list, 
#                          stats = gseagenes,
#                          minSize=3,
#                          maxSize=Inf, nperm = 100000) 
# 
# topPathwaysUp <- fgseaRes[ES>0][head(order(padj),n=15), ]
# topPathwaysDown <- fgseaRes[ES<0][head(order(padj), n=15), ]
# combipath <- rbind(topPathwaysUp, topPathwaysDown)
# combipath <- combipath %>% 
#   mutate(path4graph = str_replace(pathway, "REACTOME_", "")) %>% 
#   mutate(sens = ifelse(ES > 0, "UP", "DOWN" ))
# pathsFull_l[[k]][[CT]] <- fgseaRes
# pathsFiltered[[k]][[CT]] <- combipath

