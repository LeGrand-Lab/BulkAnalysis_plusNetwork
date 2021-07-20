# =================================================
# try to build a matrix counting the connections:
# universal colors used for cell types
# ok for color blindness!  hex codes: https://rdrr.io/cran/ggthemes/man/colorblind.html

library(circlize)
Gege = read.graph(file=paste0(grdir, "Old","_","D2",
                       "_igraph.ml"), format="graphml")

givememystuff <- function(G){
  # produces 4 matrices, for example counts:
  #         FAPs Neutro sCs ECs  M1  M2
  # FAPs    202    125 198 178 152 149
  # Neutro   57     44  67  81  36  51
  # sCs     137     71 149 112  73  79
  # ...
  tabedges = igraph::as_data_frame(G,"edges")
  tabvertx = igraph::as_data_frame(G, "vertices")
  groups = unique(igraph::vertex_attr(G, "celltype", index = V(G)))
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
    w_mat[typeorig,typedest] = w_mat[typeorig,typedest] + tabedges[i,]$weight
  }
  colormatrix =  matrix('', length(groups),length(groups) )
  rownames(colormatrix) = groups
  colnames(colormatrix) = groups
  for (i in 1:length(groups)){
    colormatrix[i,] <- rep(cellcolors[[groups[i]]], length(groups))
  }
  r_mat = w_mat/mat   # ratio matrix, weight/number
  return(list("number"=mat, "cummul_w"=w_mat, "ratio"=r_mat, "colorsmat"=colormatrix))
}

youpi = givememystuff(Gege)

# idea : make ratio button, by default circos only with RATIO, offer nb and cummul weight
cellcolors = list(
  "ECs"="#0072B2",
  "FAPs"="#F0E442",
  "M1" = "#D55E00",
  "M2" =  "#CC79A7",
  "Neutro" =  "#009E73",
  "sCs" = "#56B4E9" 
)
circos.clear()
#  add titles as thing goes plotted in apllication because here it does not work 
chordDiagram(youpi[["number"]], directional = 1, 
             grid.col = unlist(cellcolors),
             col = youpi[["colorsmat"]],
             direction.type = c("diffHeight", "arrows"),
             link.arr.type = "big.arrow",
             self.link = 1,
             annotationTrack = c("name", "grid"," axis")) 
#+ title("Number of connections")
circos.clear()

chordDiagram(youpi[["ratio"]], directional = 1, 
             grid.col = unlist(cellcolors),
             col = youpi[["colorsmat"]],
             direction.type = c("diffHeight", "arrows"),
             link.arr.type = "big.arrow",
             self.link = 1,
             annotationTrack = c("name", "grid"," axis")) 
circos.clear()
chordDiagram( youpi[["cummul_w"]], directional = 1, 
             grid.col = unlist(cellcolors),
             col = youpi[["colorsmat"]],
             direction.type = c("diffHeight", "arrows"),
             link.arr.type = "big.arrow",
             self.link = 1,
             annotationTrack = c("name", "grid"," axis")) 
circos.clear()

# circos.clear()
# chordDiagram(mat, directional = 1, 
#              grid.col = unlist(cellcolors),
#              col = colormatrix,
#              direction.type = c("diffHeight", "arrows"),
#              link.arr.type = "big.arrow",
#              self.link = 1,
#              annotationTrack = c( "grid"," axis"),
#              annotationTrackHeight = mm_h(5))
# for(si in get.all.sector.index()){
#   xlim = get.cell.meta.data("xlim", sector.index = si, track.index = 1)
#   ylim = get.cell.meta.data("ylim", sector.index = si, track.index = 1)
#   circos.text(mean(xlim),mean(ylim),si, sector.index = si,
#               track.index = 1, facing="bending.inside", niceFacing = T, col = "white")
# }

# #  ============ minimal test
# mini = matrix(c(100,200,10,50), 2,2 ) 
# rownames(mini) = c('A','B')
# colnames(mini) = c('A','B')
# mat
# #    A  B
# # A 100 10
# # ...
# circos.clear()
# chordDiagram(mini, directional = 1, 
#              grid.col = unlist(cellcolors),
#              col = colormatrix,
#              link.arr.width = 0.05,
#              link.arr.length = 0.08,
#              direction.type = "arrows",
#              link.arr.col =  "whitesmoke",           
#              self.link = 1,
#              annotationTrack = c("name", "grid"," axis")
# ) + title("my smaaaall test")
# circos.clear()




