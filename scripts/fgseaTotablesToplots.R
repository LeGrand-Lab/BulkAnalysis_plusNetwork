# Organize GSEA results into tables
# (remind: fgsea was performed on top ~100 "softfiltered" gene list
#   check fgsea_dayandcelltype.R)
# output dir: exam_INTER_conditions/static/GSEA/ 
# --
# johaGL
library(tidyverse)
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
copypathsFiltered <- pathsFiltered
for (d in names(pathsFiltered)){
  for (t in names(pathsFiltered[[d]])){
    print(dim(pathsFiltered[[d]][[t]]))
    tmp <- pathsFiltered[[d]][[t]]
    tmp <- tmp %>% mutate(pathway = sapply(path4graph, function(x){
      if (str_starts(x,"HALLMARK")){
        return(x)
      }else{return(paste0("REACTOME_",x))}
    })) %>% select(-path4graph)
    copypathsFiltered[[d]][[t]] <- tmp
  }}
saveRDS(copypathsFiltered, paste0(odir, "GSEA/fgseaByDay_filtered.rds" ) )
# TODO  : create one dfx for Reactome and another one for Hallmark 
#           and exclude DISEASE (use str_detect) and DEVELOPMENT 
patho = c()
sens = c()
NESO = c()
geneo =  c()
celltypeo = c()
dayo = c()

counter = 0
for (d in names(pathsFiltered)){
  for (t in names(pathsFiltered[[d]])){
    tmp <- pathsFiltered[[d]][[t]]
    path_ = tmp$path4graph 
    for (i in 1:length(path_)){
      pw = tmp[i,]$path4graph
      NES = tmp[i,]$NES
      se = tmp[i,]$sens
      for(lg in tmp[i,]$leadingEdge){
        if(length(lg) > 5){
          for (g in lg){
            print(g)
            counter = counter + 1
            patho <- c(patho, pw)
            sens <- c(sens, se )
            NESO <- c(NESO, NES)
            geneo <- c(geneo, g)
            celltypeo <- c(celltypeo, t)
            dayo <- c(dayo, d)
          }
        }
      }
    }
  }
}

dfx <- data.frame("path" = patho,
                  "sens" = sens,
            "NES" = NESO,
            "gene" = geneo ,
            "celltype" = celltypeo,
            "day" = dayo)
dfx <- dfx %>% mutate(unigene = paste0(gene,"_", celltype),
                      unipath = paste0(path,"_", day)) 

# build megamatrix:
mato <- array(NA,dim=c( length(unique(dfx$unipath)),
                        length(unique(dfx$unigene))
                      ))
rownames(mato) <- unique(dfx$unipath)
colnames(mato) <- unique(dfx$unigene)

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
      mato[p,g] <- val
    }
  }
}
View(mato)

keep.c <- apply(mato, 2, function(col) {
  if(length(col[is.na(col)]) >= length(col)-1){
    return(FALSE)
  }else{
    return(TRUE)
  }
} )
mati <- mato[, keep.c]
daysvec <- sapply(rownames(mati), function(x) {
  ltmp <- unlist(str_split(x, "_")) 
  return(ltmp[length(ltmp)])})
celltyvec <- sapply(colnames(mati),function(x){
  ltmp <- unlist(str_split(x, "_"))
  return(ltmp[length(ltmp)])
})

splitdays <- data.frame(x=daysvec)
rownames(splitdays) <- rownames(mati)
splitcols <- data.frame(y=celltyvec)
rownames(splitcols) <- colnames(mati)
pdf(paste0(odir,"GSEA/test1.pdf"), width = 11)
ComplexHeatmap::Heatmap(mati, na_col = "whitesmoke",
                        cluster_rows = FALSE,
                        cluster_columns = FALSE,
                        row_split = splitdays,
                        column_split = splitcols,
                        row_names_gp = gpar(fontsize = 8),
                        column_names_gp = gpar(fontsize = 8)
                        )
dev.off()


nr1 = 4; nr2 = 8; nr3 = 6; nr = nr1 + nr2 + nr3
nc1 = 6; nc2 = 8; nc3 = 10; nc = nc1 + nc2 + nc3
mati = cbind(rbind(matrix(rnorm(nr1*nc1, mean = 1,   sd = 0.5), nr = nr1),
                  matrix(rnorm(nr2*nc1, mean = 0,   sd = 0.5), nr = nr2),
                  matrix(rnorm(nr3*nc1, mean = 0,   sd = 0.5), nr = nr3)),
            rbind(matrix(rnorm(nr1*nc2, mean = 0,   sd = 0.5), nr = nr1),
                  matrix(rnorm(nr2*nc2, mean = 1,   sd = 0.5), nr = nr2),
                  matrix(rnorm(nr3*nc2, mean = 0,   sd = 0.5), nr = nr3)),
            rbind(matrix(rnorm(nr1*nc3, mean = 0.5, sd = 0.5), nr = nr1),
                  matrix(rnorm(nr2*nc3, mean = 0.5, sd = 0.5), nr = nr2),
                  matrix(rnorm(nr3*nc3, mean = 1,   sd = 0.5), nr = nr3))
)
rownames(mati) = paste0("row", seq_len(nr))
colnames(mati) = paste0("column", seq_len(nc))
library(ComplexHeatmap)
Heatmap(mati)
# --
p_ = readRDS( paste0(odir, "GSEA/fgseaByDay_filtered.rds" ) )
for (d in names(p_)){
  for (t in names(p_[[d]])){
    
  }
}

dim(dfx)

# ===================== Test some plot =============================
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