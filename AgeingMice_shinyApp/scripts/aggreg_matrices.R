# This script Bulds for each network, 4 matrices:
#  Counts, cumulatedWeight and Ratio CumulatedWeight/counts 
# which are indispensable for chord diagrams in Shiny app
# universal colors used for cell types
# ok for color blindness!  hex codes: https://rdrr.io/cran/ggthemes/man/colorblind.html
#
# johaGL 2021, Pauline
library(igraph)
<<<<<<< HEAD
library(tidyverse)

grdir <- "../graphobjs_copy/"
=======
library(stringr)

setwd("~/BulkAnalysis_plusNetwork/")
grdir <- "/networks_explore/graphobjs/"
>>>>>>> f8052c63143e8d6abf54325f61cf5382a3804a4c

dosummarymatrices <- function(G, cellcolors){
  # input :  igraph objects (plus attributes) and a compatible colors list
  # produces 4 matrices, for example counts:
  #         FAPs Neutro sCs ECs  M1  M2
  # FAPs    202    125 198 178 152 149
  # Neutro   57     44  67  81  36  51
  # sCs     137     71 149 112  73  79
  # ...
  tabedges = igraph::as_data_frame(G,"edges")
  #print(tabedges)
  tabvertx = igraph::as_data_frame(G, "vertices")
  #print(tabvertx)
  groups = unique(igraph::vertex_attr(G, "celltype", index = V(G)))
  #print(groups)
  # mat : a matrix that counts connections among given cell types 
  mat = matrix(0, length(groups),length(groups) )  # counts matrix
  rownames(mat) = groups
  colnames(mat) = groups
  w_mat = mat    #  a matrix that sums weights of those connections
  for (i in 1:dim(tabedges)[1]){
    typeorig <- unlist(str_split(tabedges[i,]$origtype, "_"))[2]
    dest <- tabedges[i,]$to
    typedest <- tabvertx[dest,]$celltype
    mat[typeorig,typedest] = mat[typeorig,typedest] + 1
    w_mat[typeorig,typedest] = round(w_mat[typeorig,typedest] + tabedges[i,]$weight, 4)
  }
  colormatrix =  matrix('', length(groups),length(groups) )
  rownames(colormatrix) = groups
  colnames(colormatrix) = groups
  for (i in 1:length(groups)){
    colormatrix[i,] <- rep(cellcolors[[groups[i]]], length(groups))
  }
  colormatrix
  r_mat = round(w_mat/mat,4)   # ratio matrix, weight/number
  return(list("Count"=mat, 
              "CumulateWeight"=w_mat, 
              "Ratio"=r_mat, 
              "Colorsmat"=colormatrix))
}
# ==================  End function ============================
# Gege = read.graph(file=paste0(grdir, "Old","_","D2",
#                      "_igraph.ml"), format="graphml")
# youpi = givememystuff(Gege)
cellcolors = list(
<<<<<<< HEAD
  "ECs"="#0072B2",
  "FAPs"="#F0E442",
  "Inflammatory-Mac" = "#D55E00",
  "Resolving-Mac" =  "#CC79A7",
  "Neutrophils" =  "#009E73",
  "MuSCs" = "#56B4E9" )  
=======
  "ECs"="#10b387ff",
  "FAPs"="#3d85c6ff",
  "MuSCs" = "#b171f1ff",
  "Neutrophils" = "#f0e442ff",
  "Inflammatory-Mac" = "#ff9900ff",
  "Resolving-Mac" = "#cc0000ff")

>>>>>>> f8052c63143e8d6abf54325f61cf5382a3804a4c
aggreg_matrices <- list()
for ( age in c("Young", "Old")){
  aggreg_matrices[[age]] <- list()
  for (day in c("D0","D2", "D4", "D7")){
<<<<<<< HEAD
    G = igraph::read.graph(file=paste0(grdir, age,"_",day, "_igraph_unfi.ml"),
                           format="graphml")
=======
    G = read.graph(file=paste0("/home/bioinfo/BulkAnalysis_plusNetwork/networks_explore/graphobjs/", age,"_",day, "_igraph_unfi.ml"), format="graphml")
>>>>>>> f8052c63143e8d6abf54325f61cf5382a3804a4c
    aggreg_matrices[[age]][[day]] <- dosummarymatrices(G, cellcolors)
  }
} 
#path="/home/bioinfo/BulkAnalysis_plusNetwork/networks_explore/graphobjs/Old_D4_igraph_unfi.ml"
#G = read.graph(file=path, format="graphml")
#aggreg_matrices[["Old"]] <- list()
#aggreg_matrices[["Old"]][["D4"]] <- dosummarymatrices(G, cellcolors)
#system("mkdir Data")
<<<<<<< HEAD
saveRDS(aggreg_matrices, "../Data/aggreg_matrices.rds")
=======
saveRDS(aggreg_matrices, "AgeingMice_shinyApp/Data/aggreg_matrices_TPM.rds")

aggreg_matrices <- list()
for ( age in c("Young", "Old")){
  aggreg_matrices[[age]] <- list()
  for (day in c("D0","D2", "D4", "D7")){
    G = read.graph(file=paste0("/home/bioinfo/BulkAnalysis_plusNetwork/networks_explore/graphobjs/", age,"_",day, "_igraph_unfi_CN.ml"), format="graphml")
    aggreg_matrices[[age]][[day]] <- dosummarymatrices(G, cellcolors)
  }
} 
saveRDS(aggreg_matrices, "AgeingMice_shinyApp/Data/aggreg_matrices_CN.rds")

>>>>>>> f8052c63143e8d6abf54325f61cf5382a3804a4c
# idea : make ratio button, by default circos only with RATIO, offer nb and cummul weight

