library(shiny)
library(igraph)
library(visNetwork)
setwd("~/BulkAnalysis_plusNetwork/graphmodeling")
##
# ATTENTION:
# only numeric and string type attributes can be retrieved
# (i.e if List or Tuples assigned to attributes from python
# they wont be read, so DON'T add that in python networkx obj)
##
g <- read_graph("myminigraph.ml", format="graphml")
# check edge attributes
(names(edge_attr(g, index=E(g))))
# origtype ecolor weight
# check node attributes
(names(vertex_attr(g,index=V(g))))
#[1] "_nx_name"  "nodetype"  "celltype"  "genesym"  
#[5] "color"     "averagexp" "id"    
is_directed(g) # TRUE

vertex_attr(g)$numid <- vertex_attr(g)$id
vertex_attr(g)$id <- vertex_attr(g)$`_nx_name`
edge_attr(g)$color <- edge_attr(g)$ecolor


l = layo
plot(g, vertex.label=V(g)$genesym,
     vertex.label.color=V(g)$color, edge.curved=.2)


#plot(g,edge.width=E(g)$weight*20,
#     sizes=(log10(V(g)$averagexp))*10, )
#layout_as_tree(g, root = numeric(), circular = FALSE,
#               rootlevel = numeric(), mode = c("all"), 
#               flip.y = TRUE)



