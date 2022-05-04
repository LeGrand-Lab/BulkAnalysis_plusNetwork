
allcombinaisonByType<-list()
sourceNode<-c()
targetNode<-c()
valueEdge<-c()
Nbconnection<-c()
for (typeCell in unique(GSEAsigni$type)){
  #keep table with signi gene
  tempoSigniDEpathway<- GSEAsigni %>% filter(type==typeCell)
  #keep only uniqueID sort by days
  uniqueIdByDay <- lapply(unique(tempoSigniDEpathway$day), function(x) filter(tempoSigniDEpathway, day == x) %>% dplyr::select(pathway) %>% unlist() %>% as.character())
  names(uniqueIdByDay) = unique(tempoSigniDEpathway$day)
  #for (d in days){ if (! d %in% names(uniqueIdByDay) ) {uniqueIdByDay[[d]]<-0 }}
  #found all combinaison where one same gene is on 1 or 2 or 3 or 4 days 
  allcombinaison =ComplexHeatmap::make_comb_mat(uniqueIdByDay)
  Sup2Combi<-allcombinaison[comb_degree(allcombinaison) >= 1]
  Node_Target<-names(comb_size(Sup2Combi))
  num_combi<-1
  for ( name in Node_Target){
    print(name)
    particulName<-str_extract_all(name,boundary("character"))[[1]]
    names(particulName)<-1:length(particulName)
    tempoName<-lapply(1:length(particulName), function(i) if(particulName[i] == "1"){daysv[[typeCell]][i]} ) %>% unlist()
    if (str_count(name,"1")==1){    Nbconnection<-c(Nbconnection,str_c(tempoName, collapse = "."))} else {Nbconnection<-c(Nbconnection,rep(Nbconnection<-c(str_c(tempoName, collapse = ".")),str_count(name,"1")-1))}
    dayDE<-c()
    countNode<-0
    firstOccurence<-T
    for (i in 1:length(particulName)){
      if (particulName[i] == "1"){
        if ( firstOccurence == T){
          sourceNode<-c(sourceNode,paste0(daysv[[typeCell]][i],"_",typeCell))
          countNode<-countNode+1
          firstOccurence<-F
          if (str_count(name,"1")==1 ){
            targetNode<-c(targetNode,paste0(daysv[[typeCell]][i],"_",typeCell))
            valueEdge<-c(valueEdge,unname(comb_size(Sup2Combi)[num_combi]))
          }
        }
        else { if (countNode<str_count(name,"1")-1){
          countNode<-countNode+1
          sourceNode<-c(sourceNode,paste0(daysv[[typeCell]][i],"_",typeCell))
          targetNode<-c(targetNode,paste0(daysv[[typeCell]][i],"_",typeCell))
          valueEdge<-c(valueEdge,unname(comb_size(Sup2Combi)[num_combi]))
        } else {
          countNode<-countNode+1
          targetNode<-c(targetNode,paste0(daysv[[typeCell]][i],"_",typeCell))
          valueEdge<-c(valueEdge,unname(comb_size(Sup2Combi)[num_combi]))
        }
        }
      }
    }
    num_combi=num_combi+ 1
  } 
}  

sourceNode<-c(sourceNode,"D0_Neutrophils","D0_Resolving-Mac2","D0_Inflammatory-Mac")
targetNode<-c(targetNode,"D2_Neutrophils","D2_Resolving-Mac","D2_Inflammatory-Mac")
valueEdge<-c(valueEdge,0,0,0)
Nbconnection<-c(Nbconnection,"D0","D0","D0")
links <- data.frame(
  source=sourceNode,
  target=targetNode,
  value=valueEdge
)
nodes <- data.frame(
  name=as.factor(c(as.character(links$source), as.character(links$target)) %>% 
                   unique())
)
nodes$group <- as.factor(c(as.character(links$source), as.character(links$target)) %>% 
                           unique())
# With networkD3, connection must be provided using id, not using real name like in the links dataframe.. So we need to reformat it.
links$IDsource <- match(links$source, nodes$name)-1 
links$IDtarget <- match(links$target, nodes$name)-1 
links$group <- as.factor(Nbconnection)
# prepare color scale: I give one specific color for each node.

my_color <- 'd3.scaleOrdinal().domain(["D0","D2","D4","D7"
                           ,"D0.D2","D2.D4","D4.D7",
                           "D0.D2.D4","D2.D4.D7","D0.D2.D4.D7",
"D0_ECs", "D2_ECs", "D0_FAPs", "D2_FAPs", "D4_FAPs","D0_MuSCs", "D2_MuSCs","D4_MuSCs","D7_MuSCs","D0_Inflammatory-Mac", "D2_Inflammatory-Mac","D4_Inflammatory-Mac","D0_Resolving-Mac" ,"D2_Resolving-Mac", "D4_Resolving-Mac","D0_Neutrophils","D2_Neutrophils"]) .range([
"#fef6be","#e2f4df","#dbeaf5","#d8d9ea",
                    "#d6efa6","#a0dbb5","#8895c2",
                    "#008349","#7e2175","#67081d",
#67001F","#67001F","#67001F","#67001F","#67001F","#67001F","#67001F","#67001F","#67001F","#67001F","#67001F","#67001F","#67001F","#67001F","#67001F","#67001F"])'

linksToRemove<-c("D0.D4","D0.D7","D2.D7","D0.D2.D4","D0.D2.D7","D0.D4.D7")
rowtoRemove<-subset(links, group %in% linksToRemove) %>% rownames()%>% as.numeric()
links2=links[-rowtoRemove,]

# Make the Network. I call my colour scale with the colourScale argument
p <- sankeyNetwork(Links = links2, Nodes = nodes, Source = "IDsource", Target = "IDtarget", 
                   Value = "value", NodeID = "name",  LinkGroup="group", NodeGroup="group", sinksRight = FALSE,
                   nodePadding=85, nodeWidth = 5,fontSize=13,fontFamily="Arial",iterations = 55,margin = list("right"=0,"left"=0), height = 642, width = 680)
p

