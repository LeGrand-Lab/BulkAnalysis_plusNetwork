# This script Bulds for each network, 4 matrices:
#  Counts, cumulatedWeight and Ratio CumulatedWeight/counts 
# which are indispensable for chord diagrams in Shiny app
# universal colors used for cell types
# ok for color blindness!  hex codes: https://rdrr.io/cran/ggthemes/man/colorblind.html
#
# johaGL 2021

setwd("~/BulkAnalysis_plusNetwork/")
grdir <- "/networks_explore/graphobjs/"

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
  "ECs"="#0072B2",
  "FAPs"="#F0E442",
  "M1" = "#D55E00",
  "M2" =  "#CC79A7",
  "Neutro" =  "#009E73",
  "sCs" = "#56B4E9" )

aggreg_matrices <- list()
for ( age in c("Young", "Old")){
  aggreg_matrices[[age]] <- list()
  for (day in c("D0","D2", "D4", "D7")){
    G = read.graph(file=paste0("/home/bioinfo/BulkAnalysis_plusNetwork/networks_explore/graphobjs/", age,"_",day, "_igraph_unfi.ml"), format="graphml")
    aggreg_matrices[[age]][[day]] <- dosummarymatrices(G, cellcolors)
  }
} 
#path="/home/bioinfo/BulkAnalysis_plusNetwork/networks_explore/graphobjs/Old_D4_igraph_unfi.ml"
#G = read.graph(file=path, format="graphml")
#aggreg_matrices[["Old"]] <- list()
#aggreg_matrices[["Old"]][["D4"]] <- dosummarymatrices(G, cellcolors)
#system("mkdir Data")
saveRDS(aggreg_matrices, "Data/aggreg_matrices.rds")
# idea : make ratio button, by default circos only with RATIO, offer nb and cummul weight

