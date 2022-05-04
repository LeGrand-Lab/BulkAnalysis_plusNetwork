
### Sankey plot by type cell

tableSigni<-signiDEgene
vectTypeCell<-"MuSCs"
outputprefix<-"MuSCs"
dblue<-"#b171f1ff"

vectTypeCell<-"ECs"
outputprefix<-"ECs"
dblue<-"#10b387ff"

vectTypeCell<-"FAPs"
outputprefix<-"FAPs"
dblue<-"#3d85c6ff"

vectTypeCell<-"Resolving-Mac"
outputprefix<-"Resolving-Mac"
dblue<-"#cc0000ff"

vectTypeCell<-"Inflammatory-Mac"
outputprefix<-"Inflammatory-Mac"
dblue<-"#ff9900ff"

### Node UP - UNDIFF - DOWN
# Nbconnection = group by Number day DEG
# NbUpDown = group by number day DEG UP and DOWN
# CombiDayDEG = group by combinaison day DEG

sourceNode<-c()
targetNode<-c()
valueEdge<-c()
NbUpDown<-c()
colorNbUpDown<-c("plum","plum4",
                 "mediumpurple1","mediumpurple4",
                 "violetred1","violetred4",
                 "LightCyan","skyblue","skyblue2","royalblue3",
                 "LightPink","indianred1","firebrick1","red3" )
colorNbUpDown2<-c("#D8BFD8","#DA70D6",
                 "#9370DB","#663399",
                 "#FF69B4","#C71585",
                 "#E0FFFF","#B0E0E6","#6495ED","#00008B",
                 "#FFB6C1","#F08080","#DC143C","#8B0000" )
names(colorNbUpDown)<-c("1_UP.1_DOWN" ,"2_UP.2_DOWN",
                        "1_UP.2_DOWN", "1_UP.3_DOWN" ,
                        "2_UP.1_DOWN","3_UP.1_DOWN",
                        "1_DOWN","2_DOWN","3_DOWN" ,"4_DOWN",
                        "1_UP.","2_UP.","3_UP.","4_UP.")
unique(sort(NbUpDown))
Nbconnection<-c()
CombiDayDEG<-c()

for ( typeCell in vectTypeCell ) {
  print(typeCell)
  #prepare list gene_SENS DEG by day
  tempoSigniDEgene<- tableSigni %>% filter(type==typeCell)
  uniqueIdByDayUP = lapply(unique(tempoSigniDEgene$day), function(x) filter(tempoSigniDEgene, day == x)  %>% filter(log2FoldChange > 0) %>% dplyr::select(id) %>% unlist() %>% as.character())
  names(uniqueIdByDayUP) = unique(paste0(tempoSigniDEgene$day,"_UP"))
  uniqueIdByDayDOWN = lapply(unique(tempoSigniDEgene$day), function(x) filter(tempoSigniDEgene, day == x)  %>% filter(log2FoldChange < 0) %>% dplyr::select(id) %>% unlist() %>% as.character())
  names(uniqueIdByDayDOWN) = unique(paste0(tempoSigniDEgene$day,"_DOWN"))
  uniqueIdByDay<-c(uniqueIdByDayUP,uniqueIdByDayDOWN)
  
  #Get combinaison
  allcombinaison =ComplexHeatmap::make_comb_mat(uniqueIdByDay)
  Sup2Combi<-allcombinaison[comb_degree(allcombinaison) >= 2]
  
  #Get name combinaison
  Node_Target<-names(comb_size(Sup2Combi))
  if(length(Node_Target) <= 2){next}
  num_combi<-1
  
  #For one combination
  #Extract the number of day DEG
  #Found the name of combination with day+sens(UP,DOWN)
  #Extract Day + sens (UP, DOWN,UNDIFF)
  for ( name in Node_Target){
    print(name)
    particulName<-str_extract_all(name,boundary("character"))[[1]]
    names(particulName)<-1:length(particulName)
    tempoName1<-lapply(1:(length(particulName)/2), function(i) if(particulName[[i]] == "1"){paste0(daysv[[typeCell]][i],"_UP")} ) %>% unlist()
    tempoName2<-lapply(((length(particulName)/2)+1):length(particulName), function(i) if(particulName[[i]] == "1"){paste0(daysv[[typeCell]][i-(length(particulName)/2)],"_DOWN")} ) %>% unlist()
    tempoName3<-lapply(1:(length(particulName)/2), function(i) if(particulName[[i]] == "0" && particulName[[i+(length(particulName)/2)]] == "0"){paste0(daysv[[typeCell]][i],"_UNDIFF")} ) %>% unlist()
    
    tempoName1bis<-lapply(1:(length(particulName)/2), function(i) if(particulName[[i]] == "1"){daysv[[typeCell]][i]} ) %>% unlist()
    tempoName2bis<-lapply(((length(particulName)/2)+1):length(particulName), function(i) if(particulName[[i]] == "1"){daysv[[typeCell]][i-(length(particulName)/2)]} ) %>% unlist()
    tempoName<-str_c(sort(c(tempoName1bis,tempoName2bis)),collapse = ".")
    
    Nodes<-sort(c(tempoName1,tempoName2,tempoName3))
    CombiDayDEG<-c(CombiDayDEG,rep(tempoName,length(Nodes)-1))
    Nbconnection<-c(Nbconnection,rep(length(c(tempoName1,tempoName2)),length(Nodes)-1))
    NbUpDown<-c(NbUpDown,rep(paste0(ifelse(length(tempoName1)==0,"",paste0(length(tempoName1),"_UP.")),ifelse(length(tempoName2)==0,"",paste0(length(tempoName2),"_DOWN"))),length(Nodes)-1))
    
    #Fill source node, target node (day consecutive) + nb gene with this combination
    sourceNode<-c(sourceNode,Nodes[1])
    if(length(Nodes)-1 > 1){
      for (i in 2:(length(Nodes)-1)){
        sourceNode<-c(sourceNode,Nodes[i])
        targetNode<-c(targetNode,Nodes[i])
        valueEdge<-c(valueEdge,unname(comb_size(Sup2Combi)[num_combi]))
      }}
    targetNode<-c(targetNode,Nodes[length(Nodes)])
    valueEdge<-c(valueEdge,unname(comb_size(Sup2Combi)[num_combi]))
    
    num_combi=num_combi+ 1
  }
}

unique(sort(NbUpDown))
links <- data.frame(
  from=str_replace_all(sourceNode,"-",""),
  to=str_replace_all(targetNode,"-",""),
  substance= as.factor(NbUpDown),
  quantity=valueEdge
)
links<-links %>% arrange(substance)
# From these flows we need to create a node data frame: it lists every entities involved in the flow
nodes <- data.frame(
  ID=c(as.character(links$from), as.character(links$to)) %>% 
    unique() %>% str_replace_all("-","")
)
nodes$label <- as.factor(c(as.character(links$from), as.character(links$to)) %>% 
                           unique())
nodes<-nodes %>% arrange(ID)
#nodes$label_pos<-c(rep("left",3),rep("above",3),rep("below",3),rep("right",3))
nodes$label_pos<-c(rep("left",3),rep("below",3),rep("right",3))
nodes$label_align<-rep("",length(nodes$ID))
nodes<-nodes %>% arrange(ID)
#nodes$x<-c(-7,"D0_DOWN+1","D0_DOWN","D0_DOWN+5","D0_DOWN+6","D0_DOWN+5","D2_DOWN+5","D2_DOWN+6","D2_DOWN+5","D4_DOWN+6","D4_DOWN+6","D4_DOWN+5")
nodes$x<-c(-7,"D2_DOWN+1","D2_DOWN","D2_DOWN+5","D2_DOWN+6","D2_DOWN+5","D4_DOWN+6","D4_DOWN+6","D4_DOWN+5")

#nodes$y<-c(-4,0,4,-5,0,5,-6,0,6,-7,0,7)
nodes$y<-c(-4,0,4,-5,0,5,-6,0,6)
nodes$dir<-rep("right",length(nodes$ID))

palettes<-data.frame(substance=names(colorNbUpDown),
                     color=colorNbUpDown)
my_title <- "Amount of genes  depending on the combination of days on a cell type \n where they are significatif differentially expressed "
attr(my_title, "gp") <- grid::gpar(fontsize=6, fontface="bold", col="black")

# node style
ns <- list(type="arrow",gp=gpar(fill=dblue, col="white", lwd=3),
           length=0.7,
           label_gp=gpar(col=dblue, fontsize=6),
           mag_pos="label", mag_fmt="%.0f", mag_gp=gpar(fontsize=6,fontface="bold",col=dblue))

pdf(paste0(analyse,"test_group_NbUpDown_NODIFF_",outputprefix,"sankeyColorUPDOWN.pdf"), width=10, height=7) # Set up PDF device
sankey(nodes, links, palettes,
       max_width=0.15, rmin=0.5,
       node_style=ns,
       page_margin=c(0, 0.05, 0, 0.05),
       title=my_title, legend=gpar(fontsize=4, col="blue", ncols=1))
dev.off()


### Sankey network
links <- data.frame(
  source=sourceNode,
  target=targetNode,
  value=sqrt(valueEdge),
  group=NbUpDown
)
links<-  aggregate(links$value,by=list("source"=links$source,"target"=links$target,"group"=links$group),FUN=sum)
# From these flows we need to create a node data frame: it lists every entities involved in the flow
nodes <- data.frame(
  name=c(as.character(links$source), as.character(links$target)) %>% 
    unique()
)
nodes$group <- as.factor(c(as.character(links$source), as.character(links$target)) %>% 
                           unique())
# With networkD3, connection must be provided using id, not using real name like in the links dataframe.. So we need to reformat it.
links$IDsource <- match(links$source, nodes$name)-1 
links$IDtarget <- match(links$target, nodes$name)-1
links$group <- as.factor(links$group)
# prepare color scale: I give one specific color for each node.
my_color <- 'd3.scaleOrdinal() .domain(["1_UP.1_DOWN","2_UP.2_DOWN","1_UP.2_DOWN","1_UP.3_DOWN","2_UP.1_DOWN","3_UP.1_DOWN","2_DOWN","3_DOWN","4_DOWN","2_UP.","3_UP.","4_UP.", 
    "D0_UP", "D2_UP","D4_UP","D7_UP", "D0_DOWN", "D2_DOWN", "D4_DOWN","D7_DOWN","D0_UNDIFF", "D2_UNDIFF","D4_UNDIFF","D7_UNDIFF"]) .range(
    ["#D8BFD8","#DA70D6",
                 "#9370DB","#663399",
                 "#FF69B4","#C71585",
                 "#B0E0E6","#6495ED","#00008B",
                 "#F08080","#DC143C","#8B0000","red","red","red","red","blue","blue","blue","blue","pink","pink","pink","pink"])'

p <- sankeyNetwork(Links = links, Nodes = nodes, Source = "IDsource", Target = "IDtarget", 
                   Value = "x", NodeID = "name", colourScale=my_color, LinkGroup="group", NodeGroup="group", sinksRight = FALSE,
                   nodePadding=50, nodeWidth = 5,fontSize=10,iterations = 55,margin = list("right"=250), height = 650, width = 850)
p
saveWidget(p, file=paste0( analyse,"test_group_NbUpDown_NODIFF_sqrt",outputprefix,"sankeyColorUPDOWN.html"))
webshot(paste0(analyse,outputprefix, "sankeyColorUPDOWN.html"), paste0(analyse,outputprefix,"sankeyCollorUPDOWN.png"), delay = 0.2)


# Without UNDIFF
# NbUpDown = group by number day DEG UP and DOWN
sourceNode<-c()
targetNode<-c()
valueEdge<-c()
NbUpDown<-c()
Nbconnection<-c()
CombiDayDEG<-c()

for ( typeCell in vectTypeCell ) {
  print(typeCell)
  #prepare list gene_SENS DEG by day
  tempoSigniDEgene<- tableSigni %>% filter(type==typeCell)
  uniqueIdByDayUP = lapply(unique(tempoSigniDEgene$day), function(x) filter(tempoSigniDEgene, day == x)  %>% filter(log2FoldChange > 0) %>% dplyr::select(id) %>% unlist() %>% as.character())
  names(uniqueIdByDayUP) = unique(paste0(tempoSigniDEgene$day,"_UP"))
  uniqueIdByDayDOWN = lapply(unique(tempoSigniDEgene$day), function(x) filter(tempoSigniDEgene, day == x)  %>% filter(log2FoldChange < 0) %>% dplyr::select(id) %>% unlist() %>% as.character())
  names(uniqueIdByDayDOWN) = unique(paste0(tempoSigniDEgene$day,"_DOWN"))
  uniqueIdByDay<-c(uniqueIdByDayUP,uniqueIdByDayDOWN)
  
  #Get combinaison
  allcombinaison =ComplexHeatmap::make_comb_mat(uniqueIdByDay)
  Sup2Combi<-allcombinaison[comb_degree(allcombinaison) >= 2]
  
  #Get name combinaison
  Node_Target<-names(comb_size(Sup2Combi))
  if(length(Node_Target) <= 2){next}
  num_combi<-1
  
  #For one combination
  #Extract the number of day DEG
  #Found the name of combination with day+sens(UP,DOWN)
  #Extract Day + sens (UP, DOWN,UNDIFF)
  for ( name in Node_Target){
    print(name)
    particulName<-str_extract_all(name,boundary("character"))[[1]]
    names(particulName)<-1:length(particulName)
    tempoName1<-lapply(1:(length(particulName)/2), function(i) if(particulName[[i]] == "1"){paste0(daysv[[typeCell]][i],"_UP")} ) %>% unlist()
    tempoName2<-lapply(((length(particulName)/2)+1):length(particulName), function(i) if(particulName[[i]] == "1"){paste0(daysv[[typeCell]][i-(length(particulName)/2)],"_DOWN")} ) %>% unlist()
    tempoName1bis<-lapply(1:(length(particulName)/2), function(i) if(particulName[[i]] == "1"){daysv[[typeCell]][i]} ) %>% unlist()
    tempoName2bis<-lapply(((length(particulName)/2)+1):length(particulName), function(i) if(particulName[[i]] == "1"){daysv[[typeCell]][i-(length(particulName)/2)]} ) %>% unlist()
    tempoName<-str_c(sort(c(tempoName1bis,tempoName2bis)),collapse = ".")
    
    Nodes<-sort(c(tempoName1,tempoName2))
    CombiDayDEG<-c(CombiDayDEG,rep(tempoName,length(Nodes)-1))
    Nbconnection<-c(Nbconnection,rep(length(c(tempoName1,tempoName2)),length(Nodes)-1))
    NbUpDown<-c(NbUpDown,rep(paste0(ifelse(length(tempoName1)==0,"",paste0(length(tempoName1),"_UP.")),ifelse(length(tempoName2)==0,"",paste0(length(tempoName2),"_DOWN"))),length(Nodes)-1))
    
    #Fill source node, target node (day consecutive) + nb gene with this combination
    sourceNode<-c(sourceNode,Nodes[1])
    if(length(Nodes)-1 > 1){
      for (i in 2:(length(Nodes)-1)){
        sourceNode<-c(sourceNode,Nodes[i])
        targetNode<-c(targetNode,Nodes[i])
        valueEdge<-c(valueEdge,unname(comb_size(Sup2Combi)[num_combi]))
      }}
    targetNode<-c(targetNode,Nodes[length(Nodes)])
    valueEdge<-c(valueEdge,unname(comb_size(Sup2Combi)[num_combi]))
    
    num_combi=num_combi+ 1
  }
}

unique(sort(NbUpDown))
links <- data.frame(
  from=str_replace_all(sourceNode,"-",""),
  to=str_replace_all(targetNode,"-",""),
  substance= as.factor(NbUpDown),
  quantity=sqrt(valueEdge)
)
links<-links %>% arrange(substance)
# From these flows we need to create a node data frame: it lists every entities involved in the flow
nodes <- data.frame(
  ID=c(as.character(links$from), as.character(links$to)) %>% 
    unique() %>% str_replace_all("-","")
)
nodes$label <- as.factor(c(as.character(links$from), as.character(links$to)) %>% 
                           unique())
nodes<-nodes %>% arrange(ID)
nodes$label_pos<-c(rep("left",2),rep("above",2),rep("below",2),rep("right",2))
nodes$label_align<-rep("",length(nodes$ID))
nodes<-nodes %>% arrange(ID)
nodes$x<-c(-7,"D0_DOWN","D0_DOWN+5","D0_DOWN+5","D2_DOWN+5","D2_DOWN+5","D4_DOWN+5","D4_DOWN+5")
nodes$y<-c(-4,4,-5,5,-6,6,-7,7)
nodes$dir<-rep("right",length(nodes$ID))

palettes<-data.frame(substance=names(colorNbUpDown),
                     color=colorNbUpDown)
my_title <- "Amount of genes (sqrt) depending on the combination of days on a cell type \n where they are significatif differentially expressed "
attr(my_title, "gp") <- grid::gpar(fontsize=6, fontface="bold", col="black")

# node style
ns <- list(type="arrow",gp=gpar(fill=dblue, col="white", lwd=3),
           length=0.7,
           label_gp=gpar(col=dblue, fontsize=6),
           mag_pos="label", mag_fmt="%.0f", mag_gp=gpar(fontsize=6,fontface="bold",col=dblue))

pdf(paste0(analyse,"test_group_NbUpDown_",outputprefix,"sankeyColorUPDOWN.pdf"), width=10, height=7) # Set up PDF device
sankey(nodes, links, palettes,
       max_width=0.15, rmin=0.5,
       node_style=ns,
       page_margin=c(0, 0.05, 0, 0.05),
       title=my_title, legend=gpar(fontsize=4, col="blue", ncols=1))
dev.off()

links <- data.frame(
  source=sourceNode,
  target=targetNode,
  value=valueEdge,
  group=Nbconnection
)
links<-  aggregate(links$value,by=list("source"=links$source,"target"=links$target,"group"=links$group),FUN=sum)
# From these flows we need to create a node data frame: it lists every entities involved in the flow
nodes <- data.frame(
  name=c(as.character(links$source), as.character(links$target)) %>% 
    unique()
)
nodes$group <- as.factor(c(as.character(links$source), as.character(links$target)) %>% 
                           unique())
# With networkD3, connection must be provided using id, not using real name like in the links dataframe.. So we need to reformat it.
links$IDsource <- match(links$source, nodes$name)-1 
links$IDtarget <- match(links$target, nodes$name)-1
links$group <- as.factor(links$group)
# prepare color scale: I give one specific color for each node.
my_color <- 'd3.scaleOrdinal() .domain(["1","2","3","4","0",
    "D0_UP", "D2_UP","D4_UP","D7_UP", "D0_DOWN", "D2_DOWN", "D4_DOWN","D7_DOWN","D0_UNDIFF", "D2_UNDIFF","D4_UNDIFF","D7_UNDIFF",]) .range(
    ["#FABDA6","#D3DDDC","#C7B19C","#446455","white","red","red","red","red","blue","blue","blue","blue","pink","pink","pink","pink"])'

}else{
  links <- data.frame(
    source=sourceNode,
    target=targetNode,
    value=valueEdge,
    group=groupetypeCell
  )
  links<-  aggregate(links$value,by=list("source"=links$source,"target"=links$target,"group"=links$group),FUN=sum)
  # From these flows we need to create a node data frame: it lists every entities involved in the flow
  nodes <- data.frame(
    name=c(as.character(links$source), as.character(links$target)) %>% 
      unique()
  )
  nodes$group <- as.factor(c(as.character(links$source), as.character(links$target)) %>% 
                             unique())
  # With networkD3, connection must be provided using id, not using real name like in the links dataframe.. So we need to reformat it.
  links$IDsource <- match(links$source, nodes$name)-1 
  links$IDtarget <- match(links$target, nodes$name)-1
  links$group <- as.factor(links$group)
  # prepare color scale: I give one specific color for each node.
  
  my_color <- 'd3.scaleOrdinal() .domain([
  "D0_UP", "D2_UP","D4_UP","D7_UP", "D0_DOWN", "D2_DOWN", "D4_DOWN","D7_DOWN","D0_UNDIFF", "D2_UNDIFF","D4_UNDIFF","D7_UNDIFF",
  "ECs","FAPs","MuSCs","Inflammatory-Mac","Resolving-Mac"]) 
  .range(["red","red","red","red","blue","blue","blue","blue","pink","pink","pink","pink","#10b387ff","#3d85c6ff","#b171f1ff","#ff9900ff","#cc0000ff"])'
  # 
  # Make the Network. I call my colour scale with the colourScale argument
  
}

p <- sankeyNetwork(Links = links, Nodes = nodes, Source = "IDsource", Target = "IDtarget", 
                   Value = "x", NodeID = "name", colourScale=my_color, LinkGroup="group", NodeGroup="group", sinksRight = FALSE,
                   nodePadding=50, nodeWidth = 5,fontSize=10,iterations = 55,margin = list("right"=250), height = 650, width = 850)
p
saveWidget(p, file=paste0( analyse,outputprefix,"sankeyColorUPDOWN.html"))
webshot(paste0(analyse,outputprefix, "sankeyColorUPDOWN.html"), paste0(analyse,outputprefix,"sankeyCollorUPDOWN.png"), delay = 0.2)

### Sankey network
links <- data.frame(
  source=sourceNode,
  target=targetNode,
  value=sqrt(valueEdge),
  group=NbUpDown
)
links<-  aggregate(links$value,by=list("source"=links$source,"target"=links$target,"group"=links$group),FUN=sum)
# From these flows we need to create a node data frame: it lists every entities involved in the flow
nodes <- data.frame(
  name=c(as.character(links$source), as.character(links$target)) %>% 
    unique()
)
nodes$group <- as.factor(c(as.character(links$source), as.character(links$target)) %>% 
                           unique())
# With networkD3, connection must be provided using id, not using real name like in the links dataframe.. So we need to reformat it.
links$IDsource <- match(links$source, nodes$name)-1 
links$IDtarget <- match(links$target, nodes$name)-1
links$group <- as.factor(links$group)
# prepare color scale: I give one specific color for each node.
my_color <- 'd3.scaleOrdinal() .domain(["1_UP.1_DOWN","2_UP.2_DOWN","1_UP.2_DOWN","1_UP.3_DOWN","2_UP.1_DOWN","3_UP.1_DOWN","2_DOWN","3_DOWN","4_DOWN","2_UP.","3_UP.","4_UP.", 
    "D0_UP", "D2_UP","D4_UP","D7_UP", "D0_DOWN", "D2_DOWN", "D4_DOWN","D7_DOWN","D0_UNDIFF", "D2_UNDIFF","D4_UNDIFF","D7_UNDIFF"]) .range(
    ["#D8BFD8","#DA70D6","#9370DB","#663399","#FF69B4","#C71585","#B0E0E6","#6495ED","#00008B","#F08080","#DC143C","#8B0000","red","red","red","red","blue","blue","blue","blue","pink","pink","pink","pink"])'

p <- sankeyNetwork(Links = links, Nodes = nodes, Source = "IDsource", Target = "IDtarget", 
                   Value = "x", NodeID = "name", colourScale=my_color, LinkGroup="group", NodeGroup="group", sinksRight = FALSE,
                   nodePadding=50, nodeWidth = 5,fontSize=10,iterations = 55,margin = list("right"=250), height = 650, width = 850)
p
saveWidget(p, file=paste0( analyse,"test_group_NbUpDown_sqrt_",outputprefix,"sankeyColorUPDOWN.html"))
webshot(paste0(analyse,outputprefix, "sankeyColorUPDOWN.html"), paste0(analyse,outputprefix,"sankeyCollorUPDOWN.png"), delay = 0.2)

#
my_color <- 'd3.scaleOrdinal().domain(["D0","D2","D4","D7"
                           ,"D0.D2","D0.D4","D0.D7","D2.D4","D2.D7","D4.D7",
                           "D0.D2.D4","D0.D2.D7","D0.D4.D7","D2.D4.D7","D0.D2.D4.D7",
"D0_ECs", "D2_ECs","D4_ECs","D7_ECs", "D0_FAPs", "D2_FAPs", "D4_FAPs","D7_FAPs","D0_MuSCs", "D2_MuSCs","D4_MuSCs","D7_MuSCs","D0_Inflammatory-Mac", "D2_Inflammatory-Mac","D4_Inflammatory-Mac","D0_Resolving-Mac" ,"D2_Resolving-Mac", "D4_Resolving-Mac","D7_Resolving-Mac","D0_Neutrophils","D2_Neutrophils"]) .range([
"#fef6be","#e2f4df","#dbeaf5","#d8d9ea",
                    "#d6efa6","#c2d9ed","#ffcce4","#a0dbb5","#d5b7d7","#8895c2",
                    "#008349","#065fa3","#d42857","#7e2175","#67081d",
"#67001F","#67001F","#67001F","#67001F","#67001F","#67001F","#67001F","#67001F","#67001F","#67001F","#67001F","#67001F","#67001F","#67001F","#67001F","#67001F","#67001F","#67001F","#67001F","#67001F","#67001F" ])'

# Prepare legend for annex to sankeyplot
colorsConditions<-c("#fef6be","#e2f4df","#dbeaf5","#d8d9ea",
                    "#d6efa6","#c2d9ed","#ffcce4","#a0dbb5","#d5b7d7","#8895c2",
                    "#008349","#065fa3","#d42857","#7e2175","#67081d")
names(colorsConditions)<-c("D0","D2","D4","D7"
                           ,"D0.D2","D0.D4","D0.D7","D2.D4","D2.D7","D4.D7",
                           "D0.D2.D4","D0.D2.D7","D0.D4.D7","D2.D4.D7","D0.D2.D4.D7")
lgddaysDEG<-Legend(at=names(colorsConditions),legend_gp = gpar(fill = colorsConditions), title="days_DEG",labels_gp = gpar (fontsize= 5),title_gp = gpar(fontsize = 6, fontface = "bold"))
lgddaysDEG_grob=grid.grabExpr(draw(lgddaysDEG)) 
grid.arrange(lgddaysDEG_grob)
###
links <- data.frame(
  from=str_replace_all(sourceNode,"-",""),
  to=str_replace_all(targetNode,"-",""),
  Nbconnection=Nbconnection,
  quantity=valueEdge
)
#Prepare links, add link between D0 and D2 if not present to have a goodposition nodes 

sourceNode<-c(sourceNode,"D0_Neutrophils","D0_Resolving-Mac","D0_Inflammatory-Mac")
targetNode<-c(targetNode,"D2_Neutrophils","D2_Resolving-Mac","D2_Inflammatory-Mac")
valueEdge<-c(valueEdge,0,0,0)
Nbconnection<-c(Nbconnection,"D0","D0","D0")
links <- data.frame(
  source=sourceNode,
  target=targetNode,
  value=sqrt(valueEdge)
)


# From these flows we need to create a node data frame: it lists every entities involved in the flow
nodes <- data.frame(
  name=c(links$source, links$target) %>% unique()
)
nodes$group <- as.factor(c(as.character(links$source), as.character(links$target)) %>% 
                           unique())
# With networkD3, connection must be provided using id, not using real name like in the links dataframe.. So we need to reformat it.
links$IDsource <- match(links$source, nodes$name)-1 
links$IDtarget <- match(links$target, nodes$name)-1
links$group <- as.factor(Nbconnection)
linksToRemove<-c("D0.D4","D0.D7","D2.D7","D0.D2.D7","D0.D4.D7")
rowtoRemove<-subset(links, group %in% linksToRemove) %>% rownames()%>% as.numeric()
links2=links[-rowtoRemove,]
# Make the Network. I call my colour scale with the colourScale argument
p <- sankeyNetwork(Links = links2, Nodes = nodes, Source = "IDsource", Target = "IDtarget", 
                   Value = "value", NodeID = "name", colourScale=my_color, LinkGroup="group", NodeGroup="group", sinksRight = FALSE,
                   nodePadding=85, nodeWidth = 5,fontSize=13,fontFamily="Arial",iterations = 55,margin = list("right"=0,"left"=0), height = 642, width = 680)
p