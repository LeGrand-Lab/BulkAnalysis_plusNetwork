# Organize GSEA results into tables
# (remind: fgsea was performed on top ~100 "softfiltered" gene list
#   check fgsea_dayandcelltype.R)
# output dir: exam_INTER_conditions/static/GSEA/ 
# --
# johaGL
library(tidyverse)
setwd("~/BulkAnalysis_plusNetwork/")
odir <- "exam_INTER_conditions/static/GSEA/"
fullrds <- readRDS(paste0(odir,"fgseaByDay_full.rds"))
tablesneeded <- F
mat4heatmapneeded <- F


if (tablesneeded){
  system(paste("cd",odir,";", 
               "if [ ! -d csv/ ]; then mkdir csv; else echo 'csv/ exists, no action';fi"))
  print(" ** see keys in full results rds list of lists ** ")
  print(names(fullrds)) # the days
  print(names(fullrds[["D0"]])) # the cell types (that exist in given day)
  days = names(fullrds)
  for (d in names(fullrds)){
    print(paste("\n^^^^^", d, "^^^^^"))
    for (t in names(fullrds[[d]])){
      print(paste("saving", d, t, "as csv"))
      print(t)
      print(dim(fullrds[[d]][[t]]))
      tmpdf <- fullrds[[d]][[t]]
      
      tmpdf <- tmpdf %>% 
        mutate(gene_symbols = sapply(leadingEdge, function(x) paste(x, collapse=", ")),
               .before = leadingEdge) %>%
        select(-leadingEdge)
      
      write.table(tmpdf, paste0(odir,"csv/", t,"_",d,"_pathways.csv"), sep="\t", 
                  col.names = T, row.names=F              )
    }
  }
}
# ===================== Prepare matrices for plots =============================
pathsFiltered = readRDS( paste0(odir, "fgseaByDay_filtered.rds" ) )

## prepare matrices for heatmaps
if (mat4heatmapneeded){
  mhp = list()
  for (k in c('D0','D2', 'D4', 'D7')){
    cts <- unique(names(pathsFiltered[[k]])) 
    mhp[[k]] <- list()
    for (CT in cts){
      mhp[[k]][[CT]] <- list("UP"=NULL,"DOWN"=NULL)
      for(S in c("UP","DOWN")){
        here.data = as.data.frame(pathsFiltered[[k]][[CT]])
        here.data <- here.data %>% filter(sens == S)
        here.dico = list()
        if(length(here.data$path4graph) == length(unique(here.data$path4graph) )){
          print(paste("ok, processing", k, CT, S))
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
          matrx <- matrx[rownames(matrx) != "DISEASE",] # take away "DISEASE" ! 
          rtokeep <- apply(matrx, 1, function(x) sum(!is.na(x)) >= 3)
          newm <- matrx[rtokeep,]
          mhp[[k]][[CT]][[S]] <- newm
        } else {
          print(paste("error, this", k, CT ,"  contains repeated pathwaynames") )
          stop()
        } # endif pathways vectors are unique (no elements repeated)
      } # end for sens UP or DOWN
    }# end for CT in cts
  }# end for k in vector of days
  saveRDS(mhp, file=paste0(odir,"GSEA/fgsea_matrices4heatmaps.rds" ))
} # end if mat4heatmapneeded
# ================================================================================
  
# ===================== Prepare giant single matrix =============================
mhp = readRDS(paste0(odir,"GSEA/fgsea_matrices4heatmaps.rds" ))

# ================================================================================
# ===================== Test some plot =============================
mhp = readRDS(paste0(odir,"GSEA/fgsea_matrices4heatmaps.rds" ))
library(pheatmap)
library(cowplot)
library(heatmap3)
library(RColorBrewer)
library(gridExtra)

names(mhp)
names(mhp[["D7"]])
# test on single mapcolor with up and down pathways
k = "D4"
CT = "M2"
pdf(paste0(odir,"t.pdf"), width=8, height=5)
heatmap3(as.matrix(mhp[[k]][[CT]][["UP"]]),  
         scale = "none",
         col=colorRampPalette(c("gray","firebrick3"))(256), 
         Colv = NA, Rowv = NA,
         cexRow = 0.7,
         cexCol = 0.7,
         ColSideWidth = ncol(mhp[[k]][[CT]][["UP"]]),
         margins = c(10,10),
         main = paste(k, CT, "UP"))

heatmap3(as.matrix(mhp[[k]][[CT]][["DOWN"]]),  
         scale = "none",
         col=colorRampPalette(c("navy", "gray"))(256), 
         Colv = NA, Rowv = NA,
         cexRow = 0.7,
         cexCol = 0.7,
         ColSideWidth = ncol(mhp[[k]][[CT]][["DOWN"]]),
         margins = c(10,10),
         main = paste(k, CT, "DOWN"))
dev.off()
print("end")