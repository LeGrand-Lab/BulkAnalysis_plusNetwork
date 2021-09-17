# Organize GSEA results into tables
# (remind: fgsea was performed on top ~100 "softfiltered" gene list
#   check fgsea_dayandcelltype.R)
# output dir: exam_INTER_conditions/static/GSEA/ 
# --
# johaGL
library(tidyverse)
library(ComplexHeatmap)

setwd("~/BulkAnalysis_plusNetwork/")
odir <- "exam_INTER_conditions/static/"
tablesneeded <- F
mat4heatmapneeded <- F

fullrds <- readRDS(paste0(odir,"GSEA/fgseaByDay_full.rds"))
pathsFiltered = readRDS( paste0(odir, "GSEA/fgseaByDay_filtered.rds" ) )

if (tablesneeded){
  system(paste0("cd ",odir,"GSEA/; ", 
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
# ============= Prepare matrices for individual (dummy) plots ==================
## prepare matrices for heatmaps
if (mat4heatmapneeded){
  m_ = list()
  for (k in c('D0','D2', 'D4', 'D7')){
    cts <- unique(names(pathsFiltered[[k]])) 
    m_[[k]] <- list()
    for (CT in cts){
      m_[[k]][[CT]] <- list("UP"=NULL,"DOWN"=NULL)
      for(S in c("UP","DOWN")){
        here.data = as.data.frame(pathsFiltered[[k]][[CT]])
        here.data <- here.data %>% filter(sens == S) %>% 
          mutate(path4graph = str_replace(pathway, "REACTOME_", "")) 
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
          m_[[k]][[CT]][[S]] <- newm
        } else {
          print(paste("error, this", k, CT ,"  contains repeated pathwaynames") )
          stop()
        } # endif pathways vectors are unique (no elements repeated)
      } # end for sens UP or DOWN
    }# end for CT in cts
  }# end for k in vector of days
  saveRDS(m_, file=paste0(odir,"GSEA/fgsea_matrices4heatmaps.rds" ))
} # end if mat4heatmapneeded

# ===================== Prepare giant single matrix ============================

# TODO  : create one dfx for Reactome and another one for Hallmark 
#           and exclude DISEASE (use str_detect) and DEVELOPMENT 
r_path = c(); r_sens = c(); r_NES = c(); r_gene = c(); r_celltype = c(); r_day = c()
h_path = c(); h_sens = c(); h_NES = c(); h_gene = c(); h_celltype = c(); h_day = c()

# split into two dataframes : one for reactome, one for hallmark
counter = 0
cutoffnbg <- 3
for (d in names(pathsFiltered)){
  for (t in names(pathsFiltered[[d]])){
    tmp <- pathsFiltered[[d]][[t]]
    path_ = tmp$path
    for (i in 1:length(path_)){
      pw = tmp[i,]$path
      NES = tmp[i,]$NES
      se = tmp[i,]$sens
      if (str_detect(pw, "DISEASE") | str_detect(pw, "DEVELOPMENT")){
        print("excluded disease and devel row")
      } else if (str_detect(pw, "REACTOME")){
        for(lg in tmp[i,]$leadingEdge){
          if(length(lg) > cutoffnbg){
            for (g in lg){
              print(g)
              counter = counter + 1
              r_path = c(r_path, pw)
              r_sens = c(r_sens, se) 
              r_NES = c(r_NES, NES)
              r_gene =  c(r_gene, g)
              r_celltype = c(r_celltype, t)
              r_day = c(r_day, d)
            } 
          }
        }
      }else if (str_detect(pw, "HALLMARK")){
        for(lg in tmp[i,]$leadingEdge){
          if(length(lg) > cutoffnbg){
            for (g in lg){
              print(g)
              counter = counter + 1
              h_path = c(h_path, pw)
              h_sens = c(h_sens, se) 
              h_NES = c(h_NES, NES)
              h_gene =  c(h_gene, g)
              h_celltype = c(h_celltype, t)
              h_day = c(h_day, d)
            } 
          }
        }
      }
    }#end for i in 1:length(path_)
  }
}

df_reactome <- data.frame("path" = r_path, "sens"= r_sens, "NES" = r_NES,
                          "gene" = r_gene, "celltype" = r_celltype, "day"=r_day)
df_hallmark <- data.frame("path" = h_path, "sens"= h_sens, "NES" = h_NES,
                          "gene" = h_gene, "celltype" = h_celltype, "day"=h_day)

df_reactome <- df_reactome %>% mutate(unigene = paste0(gene,"_", celltype),
                      unipath = paste0(path,"_", day)) 
df_hallmark <- df_hallmark %>% mutate(unigene = paste0(gene,"_", celltype),
                                      unipath = paste0(path,"_", day)) 

# build megamatrix:

dfTomatrix <- function(dfx){
  
  mat <- array(NA,dim=c( length(unique(dfx$unipath)),
                          length(unique(dfx$unigene))
  ))
  rownames(mat) <- unique(dfx$unipath)
  colnames(mat) <- unique(dfx$unigene)
  
  for (p in dfx$unipath){
    for (g in dfx$unigene){
      tryCatch({
        val <- dfx[dfx$unipath == p & dfx$unigene == g,]$NES
      },
      warning = function(w){
      },
      error = function(e){
      })
      if (length(val)!= 0 ){
        mat[p,g] <- val
      }
    }
  }
  return(mat)
}

m_reac <- dfTomatrix(df_reactome)
m_hall <- dfTomatrix(df_hallmark)

clearmatrix <- function(mat, m){
  keep.c <- apply(mat, 2, function(col) {
    if(length(col[is.na(col)]) >= length(col)-1){
      return(FALSE) # at least 1 NES values by gene
    }else{ return(TRUE) } } )
  mati <- mat[, keep.c]
  keep.r <- apply(mati, 1, function(row){
    if(length(row[is.na(row)]) >= length(row)-m){
      return(FALSE) # at least m NES values by path
    }else{ return(TRUE) } } )
  mati <- mati[keep.r, ]
  keep.cbis <- apply(mati, 2, function(col){
    if(length(col[is.na(col)]) >= length(col)-1){
      return(FALSE)
    }else{ return(TRUE) }
  })
  mati <- mati[, keep.cbis]
  return(mati)
}

m_reac2 <- clearmatrix(m_reac, 3)
m_hall2 <- clearmatrix(m_hall, 3)
# =============================== REACTOME plot ================================== 
# render plot reactome

reac_daysvec <- sapply(rownames(m_reac2), function(x) {
  ltmp <- unlist(str_split(x, "_")) 
  return(ltmp[length(ltmp)])})
reac_celltyvec <- sapply(colnames(m_reac2),function(x){
  ltmp <- unlist(str_split(x, "_"))
  return(ltmp[length(ltmp)])
})
reac_newrownames <-  make.unique( sapply( names(reac_daysvec), function(x){  
  # "REACTOME_COLLAGEN_FORMATION_D7"
  ltmp <- unlist(str_split(str_replace(x,"REACTOME_",""), "_"))  
  mot <- paste0(ltmp[1:(length(ltmp)-1)], collapse="_")
  return(mot)  # "COLLAGEN_FORMATION
}) )
reac_newcolnames <- make.unique( sapply(names(reac_celltyvec), function(x){
  ltmp <- unlist(str_split(x, "_"))  # "Lamc1_sCs"
  return(ltmp[1]) # Lamc1
}) )
# assigne new names to all objects
rownames(m_reac2) <- reac_newrownames
colnames(m_reac2) <- reac_newcolnames
names(reac_daysvec) <- reac_newrownames
names(reac_celltyvec) <- reac_newcolnames

reac_splitdays <- data.frame(x=reac_daysvec)
rownames(reac_splitdays) <- rownames(m_reac2)
reac_splitcols <- data.frame(y=reac_celltyvec)
rownames(reac_splitcols) <- colnames(m_reac2)

ha = HeatmapAnnotation(foo = anno_block(gp = gpar(fill = c("darkblue", "orange",
                                                           "firebrick","violet",
                                                           "darkgreen","royalblue") )
)) #  labels = c("ECs","FAPs","M1","M2","Neutro","sCs") 
oh <- ComplexHeatmap::Heatmap(m_reac2, 
                              na_col = "floralwhite",
                              cluster_rows = FALSE,
                              cluster_columns = FALSE,
                              row_split = reac_splitdays,
                              column_split = reac_splitcols,
                              row_names_gp = gpar(fontsize = 8),
                              column_names_gp = gpar(fontsize = 8),
                              name = "REACTOME",
                              column_names_rot = 45,
                              heatmap_legend_param = list(title="NES",
                                                          direction = "horizontal"),
                              top_annotation = ha
)

pdf(paste0(odir,"GSEA/test4.pdf"), width = 20, height = 13)

draw (oh, column_title = "GSEA Old vs Young (REACTOME)",
      heatmap_legend_side = "bottom", padding = unit(c(2, 2, 2, 60), "mm"))

dev.off()

# ================================== HALLMARK ================================== 
# TODO :  h_

# --


# ===================== Test some other plot =============================
m_ = readRDS(paste0(odir,"GSEA/fgsea_matrices4heatmaps.rds" ))
library(pheatmap)
library(cowplot)
library(heatmap3)
library(RColorBrewer)
library(gridExtra)

names(m_)
names(m_[["D7"]])
# test on single mapcolor with up and down pathways
k = "D4"
CT = "M2"
pdf(paste0(odir,"t.pdf"), width=8, height=5)
heatmap3(as.matrix(m_[[k]][[CT]][["UP"]]),  
         scale = "none",
         col=colorRampPalette(c("gray","firebrick3"))(256), 
         Colv = NA, Rowv = NA,
         cexRow = 0.7,
         cexCol = 0.7,
         ColSideWidth = ncol(m_[[k]][[CT]][["UP"]]),
         margins = c(10,10),
         main = paste(k, CT, "UP"))

heatmap3(as.matrix(m_[[k]][[CT]][["DOWN"]]),  
         scale = "none",
         col=colorRampPalette(c("navy", "gray"))(256), 
         Colv = NA, Rowv = NA,
         cexRow = 0.7,
         cexCol = 0.7,
         ColSideWidth = ncol(m_[[k]][[CT]][["DOWN"]]),
         margins = c(10,10),
         main = paste(k, CT, "DOWN"))
dev.off()
print("end")



# for (d in names(pathsFiltered)){
#   for (t in names(pathsFiltered[[d]])){
#     tmp <- pathsFiltered[[d]][[t]]
#     path_ = tmp$path 
#     if (str_detect(path_, "DISEASE") | str_detect(path_, "DEVELOPMENT")){
#       print("excluded DISEASE or DEVELOMPENT row")
#     }else if (str_detect(path_, "REACTOME")){
#       for (i in 1:length(path_)){
#         pw = tmp[i,]$path4graph
#         NES = tmp[i,]$NES
#         se = tmp[i,]$sens
#         for(lg in tmp[i,]$leadingEdge){
#           if(length(lg) > 5){
#             for (g in lg){
#               print(g)
#               counter = counter + 1
#               r_path = c(r_path, pw)
#               r_sens = c(r_sens, se) 
#               r_NES = c(r_NES, NES)
#               r_gene =  c(r_gene, g)
#               r_celltype = c(r_celltype, t)
#               r_day = c(r_day, d)
#               
#             } 
#           }
#         }
#       } #end for i in 1:length(path_)
#     } else if (str_detect(path_, "HALLMARK")){
#       for (i in 1:length(path_)){
#         pw = tmp[i,]$path4graph
#         NES = tmp[i,]$NES
#         se = tmp[i,]$sens
#         for(lg in tmp[i,]$leadingEdge){
#           if(length(lg) > 5){
#             for (g in lg){
#               print(g)
#               counter = counter + 1
#               h_path = c(h_path, pw)
#               h_sens = c(h_sens, se) 
#               h_NES = c(h_NES, NES)
#               h_gene =  c(h_gene, g)
#               h_celltype = c(h_celltype, t)
#               h_day = c(h_day, d)
#               
#             } 
#           }
#         }
#       } #end for i in 1:length(path_)
#     }
#   }
# }


# copypathsFiltered <- pathsFiltered
# for (d in names(pathsFiltered)){
#   for (t in names(pathsFiltered[[d]])){
#     print(dim(pathsFiltered[[d]][[t]]))
#     tmp <- pathsFiltered[[d]][[t]]
#     tmp <- tmp %>% mutate(pathway = sapply(path4graph, function(x){
#       if (str_starts(x,"HALLMARK")){
#         return(x)
#       }else{return(paste0("REACTOME_",x))}
#     })) %>% select(-path4graph)
#     copypathsFiltered[[d]][[t]] <- tmp
#   }}
# saveRDS(copypathsFiltered, paste0(odir, "GSEA/fgseaByDay_filtered.rds" ) )



# ComplexHeatmap::Heatmap(m_reac2, na_col = "whitesmoke",
#                         cluster_rows = FALSE,
#                         cluster_columns = FALSE,
#                         row_split = reac_splitdays,
#                         column_split = reac_splitcols,
#                         row_names_gp = gpar(fontsize = 8),
#                         column_names_gp = gpar(fontsize = 8),
#                         name = "REACTOME",
#                         column_names_rot = 45, 
#                         heatmap_legend_param = list(title="NES", 
#                                                     just = c("left", "bottom"))
# )