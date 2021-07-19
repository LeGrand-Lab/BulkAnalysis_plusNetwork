createlayoutbygroup <- function(G, group=celltype){
  G_Grouped = G
  
  return(G_Grouped)
}
age = 'Old'

G = read.graph(file=paste0(grdir, age,"_", "D0",
                           "_igraph.ml"), format="graphml")
G_Grouped = G
E(G_Grouped)$weight = 1


V(G)$group = igraph::vertex_attr(G, "celltype", index = V(G))
groups = unique(igraph::vertex_attr(G, "celltype", index = V(G)))

E(G)$color = igraph::edge_attr(G,"ecolor", index=E(G))

for (i in groups){
  GroupV = which(V(G)$group == i)
  G_Grouped = add_edges(G_Grouped,combn(GroupV,2), attr=list(weight=3))
}

LO = layout_with_fr(G_Grouped)
plot(G,vertex.size=5, mark.border=NA,
     edge.border=NA,  vertex.border="white",
     vertex.label =NA, edge.arrow.size=0.005,
     edge.arrow.width=0.00001, layout=LO)
