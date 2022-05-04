
library(grid)
library(PantaRhei)
#tableSankeyplot function to extrat table nodes sources and target with value , for a sankey plot
#plot gathered gene with the same pattern (for one type cell, witch days are DEG and in witch sens) and links days in function of sens of log2foldchange
#so if gene is UP, it is more exessed in OLD condition than Young
#signiDEgene = table with all DEG padj<0.05
#vectTypeCell = type cell selected, vector possible
#minCombi = founds genes with at least 1 or 2 days DEG
#NODIFF = Nodes inclus days UNDIFF

tableSankeyplot<- function(signiDEgene,vectTypeCell,minCombi,NODIFF){
sourceNode<-c()
targetNode<-c()
valueEdge<-c()
NbUpDown<-c()
Nbconnection<-c()
CombiDayDEG<-c()

for ( typeCell in vectTypeCell ) {
  print(typeCell)
  #prepare list gene_SENS DEG by day
  tempoSigniDEgene<- signiDEgene %>% filter(type==typeCell)
  uniqueIdByDayUP = lapply(unique(tempoSigniDEgene$day), function(x) filter(tempoSigniDEgene, day == x)  %>% filter(log2FoldChange > 0) %>% dplyr::select(id) %>% unlist() %>% as.character())
  names(uniqueIdByDayUP) = unique(paste0(tempoSigniDEgene$day,"_UP"))
  uniqueIdByDayDOWN = lapply(unique(tempoSigniDEgene$day), function(x) filter(tempoSigniDEgene, day == x)  %>% filter(log2FoldChange < 0) %>% dplyr::select(id) %>% unlist() %>% as.character())
  names(uniqueIdByDayDOWN) = unique(paste0(tempoSigniDEgene$day,"_DOWN"))
  uniqueIdByDay<-c(uniqueIdByDayUP,uniqueIdByDayDOWN)
  
  #Get combinaison
  allcombinaison =ComplexHeatmap::make_comb_mat(uniqueIdByDay)
  Sup2Combi<-allcombinaison[comb_degree(allcombinaison) >= minCombi]
  
  #Get name combinaison
  Node_Target<-names(comb_size(Sup2Combi))
  if(length(Node_Target) < minCombi){next}
  num_combi<-1
  
  #For one combination
  #Extract the number of day DEG
  #Found the name of combination with day+sens(UP,DOWN)
  #Extract Day + sens (UP, DOWN,UNDIFF)
  for ( name in Node_Target){
    print(name)
    particulName<-str_extract_all(name,boundary("character"))[[1]]
    names(particulName)<-1:length(particulName)
    #extract day DEG UP
    tempoName1<-lapply(1:(length(particulName)/2), function(i) if(particulName[[i]] == "1"){paste0(daysv[[typeCell]][i],"_UP")} ) %>% unlist()
    #extract day DEG DOWN
    tempoName2<-lapply(((length(particulName)/2)+1):length(particulName), function(i) if(particulName[[i]] == "1"){paste0(daysv[[typeCell]][i-(length(particulName)/2)],"_DOWN")} ) %>% unlist()
    #Extract day UNDIFF
    tempoName3<-lapply(1:(length(particulName)/2), function(i) if(particulName[[i]] == "0" && particulName[[i+(length(particulName)/2)]] == "0"){paste0(daysv[[typeCell]][i],"_UNDIFF")} ) %>% unlist()
    
    #Extract only day DEG
    tempoName1bis<-lapply(1:(length(particulName)/2), function(i) if(particulName[[i]] == "1"){daysv[[typeCell]][i]} ) %>% unlist()
    tempoName2bis<-lapply(((length(particulName)/2)+1):length(particulName), function(i) if(particulName[[i]] == "1"){daysv[[typeCell]][i-(length(particulName)/2)]} ) %>% unlist()
    tempoName<-str_c(sort(c(tempoName1bis,tempoName2bis)),collapse = ".")
    
    if (NODIFF == T ){
    Nodes<-sort(c(tempoName1,tempoName2,tempoName3))
    CombiDayDEG<-c(CombiDayDEG,rep(tempoName,length(Nodes)-1))
    Nbconnection<-c(Nbconnection,rep(length(c(tempoName1,tempoName2)),length(Nodes)-1))
    NbUpDown<-c(NbUpDown,rep(paste0(ifelse(length(tempoName1)==0,"",paste0(length(tempoName1),"_UP.")),ifelse(length(tempoName2)==0,"",paste0(length(tempoName2),"_DOWN"))),length(Nodes)-1))
    
    }else{    
    Nodes<-sort(c(tempoName1,tempoName2))
    if (length(Nodes)==1){Nbconnection<-c(Nbconnection,rep(length(c(tempoName1,tempoName2)),4))} else {Nbconnection<-c(Nbconnection,rep(length(c(tempoName1,tempoName2)),length(Nodes)-1))}
    if (length(Nodes)==1){CombiDayDEG<-c(CombiDayDEG,rep(tempoName,4))} else {CombiDayDEG<-c(CombiDayDEG,rep(tempoName,length(Nodes)-1))}
    if (length(Nodes)==1){NbUpDown<-c(NbUpDown,rep(paste0(ifelse(length(tempoName1)==0,"",paste0(length(tempoName1),"_UP.")),ifelse(length(tempoName2)==0,"",paste0(length(tempoName2),"_DOWN"))),4))} else {NbUpDown<-c(NbUpDown,rep(paste0(ifelse(length(tempoName1)==0,"",paste0(length(tempoName1),"_UP.")),ifelse(length(tempoName2)==0,"",paste0(length(tempoName2),"_DOWN"))),length(Nodes)-1))}
    
    }

    #Fill source node, target node (day consecutive) + nb gene with this combination
    if (length(Nodes) == 1){
        sourceNode<-c(sourceNode,Nodes[1])
        targetNode<-c(targetNode,paste0(Nodes[1],"_only1"))
        valueEdge<-c(valueEdge,unname(comb_size(Sup2Combi)[num_combi]))
        sourceNode<-c(sourceNode,paste0(Nodes[1],"_only1"))
        targetNode<-c(targetNode,paste0(Nodes[1],"_only2"))
        valueEdge<-c(valueEdge,unname(comb_size(Sup2Combi)[num_combi]))
        sourceNode<-c(sourceNode,paste0(Nodes[1],"_only2"))
        targetNode<-c(targetNode,paste0(Nodes[1],"_only3"))
        valueEdge<-c(valueEdge,unname(comb_size(Sup2Combi)[num_combi]))
        sourceNode<-c(sourceNode,paste0(Nodes[1],"_only3"))
        targetNode<-c(targetNode,Nodes[1])
        valueEdge<-c(valueEdge,unname(comb_size(Sup2Combi)[num_combi]))
    }else{
    sourceNode<-c(sourceNode,Nodes[1])
    if(length(Nodes)-1 > 1){
      for (i in 2:(length(Nodes)-1)){
        sourceNode<-c(sourceNode,Nodes[i])
        targetNode<-c(targetNode,Nodes[i])
        valueEdge<-c(valueEdge,unname(comb_size(Sup2Combi)[num_combi]))
      }}
    targetNode<-c(targetNode,Nodes[length(Nodes)])
    valueEdge <-c(valueEdge,unname(comb_size(Sup2Combi)[num_combi]))
    }
    
    num_combi=num_combi+ 1
  }

}

links <- data.frame(
  from=str_replace_all(sourceNode,"-",""),
  to=str_replace_all(targetNode,"-",""),
  NbUpDown= NbUpDown,
  Nbconnection=Nbconnection,
  CombiDayDEG=CombiDayDEG,
  quantity=valueEdge
)


return(links)
}

# sankeyPantaRhei , made a sankey plot, where you can choice coordinate of nodes, label nodes,group color links by group (NbUpDown,Nbconnection,CombiDayDEG)
#links1 = output of tableSankeyplot with at least nodes from, nodes to, value of the flow, and group
#dblue= color of nodes
#group=
#NbUpDown: group links by pattern number of UP and number od DOWN
#Nbconnection: group by number of day DEG
#CombiDayDEG: group by pattern of days DEG
#colorvect= vector color with names unique(links1$group)
#outputprefix= start of plot title
#sqrt= T/F, value in square root
#NODIFF= T/F , nodes UNDIFF?
sankeyPantaRhei<-function(links1,dblue,group,colorvect,outputprefix,sqrt,NODIFF){
  if ( sqrt == T){
    links1$quantity<-sqrt(links1$quantity)
    titlesqrt<-"(sqrt)"}else{titlesqrt<-""
    }
  print(links1)
  links <- data.frame(
    from=links1$from,
    to=links1$to,
    substance= as.factor(links1 %>% dplyr::select(group) %>% unlist()),
    #substance= as.factor(links1$CombiDayDEG),
    quantity=links1$quantity
  )
  print(links)
  links<-links %>% arrange(substance)
  # From these flows we need to create a node data frame: it lists every entities involved in the flow
  nodes <- data.frame(
    ID=c(as.character(links$from), as.character(links$to)) %>% 
      unique() %>% str_replace_all("-","")
  )
  nodes$label <- as.factor(c(as.character(links$from), as.character(links$to)) %>% 
                             unique())
  nodes<-nodes  %>% mutate(day=sapply(strsplit(nodes$ID,"_"), `[`, 1)) %>% mutate(sens=sapply(strsplit(nodes$ID,"_"), `[`, 2)) %>% mutate(loops=sapply(strsplit(nodes$ID,"_"), `[`, 3))
  nbDay<-length(unique(nodes$day))
  print(nbDay)
  nodes<-nodes  %>% arrange(loops,sens,day)
  print(nodes)
  if("only1" %in% nodes$loops){only=T}else{only=F}
  if(length(nodes %>% filter(sens== "UNDIFF") %>% dplyr::select(sens) %>% unlist()) >=1){
    replabel=3
    titleNODIFF="UNDIFF"}else{
    NODIFF = F
    replabel=2
    titleNODIFF=""}
  if (nbDay == 4){
    poslabel=c("left","above","below","right")
    xpos<-c(-7,-2,3,8)
    ypos<-c(-2,-3,-4,-5)
  } else if(nbDay == 3){
    poslabel=c("left","above","right")
    xpos<-c(-7,-2,3)
    ypos<-c(-2,-3,-4)
  } else{
      poslabel=c("left","right")
      xpos<-c(-7,3)
      ypos<-c(-2,-3)
  }

  nodes<- nodes %>% mutate(label=lapply(1:length(nodes$sens) ,function(i) ifelse( is.na(nodes[i,]$loops) == T,nodes[i,]$ID,"" )) %>% unlist()) 
  nodes<- nodes %>%  mutate(ID=lapply( 1:length(nodes$sens) ,function(i) ifelse( is.na(nodes[i,]$loops) == F,paste0(".",nodes[i,]$ID), nodes[i,]$ID)) %>% unlist())
  if(only == T){nodes$label_pos<-c(rep("none",nbDay*6),rep(poslabel,replabel))}else{nodes$label_pos<-c(rep(poslabel,replabel))}
  nodes$label_align<-rep("",length(nodes$ID))
  if(NODIFF == T){nodes$x<-c(xpos,xpos+1,xpos)}else{ if(only == T){nodes$x<-c(xpos+2,xpos+2,xpos,xpos,xpos-2,xpos-2,rep(xpos,2))} else{nodes$x<-c(xpos,xpos)}}
  if(NODIFF == T){nodes$y<-c(ypos,rep(0,nbDay),ypos*-1)}else{ if(only == T){nodes$y<-c(ypos-1,(ypos*-1)+1,ypos-2,(ypos*-1)+2,ypos-1,(ypos*-1)+1,ypos,ypos*-1)}else{nodes$y<-c(ypos,ypos*-1)}}
  if(only == T){nodes$dir<-c(rep("down",nbDay),rep("up",nbDay),rep("left",nbDay*2),rep("up",nbDay),rep("down",nbDay),rep("right",nbDay*2))}else{nodes$dir<-rep("right",nbDay*replabel)}
  
  palettes<-data.frame(substance=names(colorvect),
                       color=colorvect)
  my_title <- paste0(outputprefix," group by ",group," with sens UP DOWN ",titleNODIFF,"\n Amount of gene ",titlesqrt," depending on the combination of days on a cell type , padj<0.005 ")
  attr(my_title, "gp") <- grid::gpar(fontsize=6, fontface="bold", col="black")
  
  # node style
  ns <- list(type="arrow",gp=gpar(fill=dblue, col="white", lwd=3),
             length=0.5,
             label_gp=gpar(col=dblue, fontsize=6,fontface="bold"),
             mag_pos="label", mag_fmt="%.0f", mag_gp=gpar(fontsize=6,fontface="bold",col=dblue))
  
  #pdf(paste0(analyse,"test_group_",group,"_",titleNODIFF,'_',sqrt,'_',outputprefix,"_sankeyColorUPDOWN.pdf"), width=10, height=7) # Set up PDF device
  print(nodes)
  PantaRhei::sankey(nodes, links, palettes,
         max_width=0.15, rmin=0.5,
         node_style=ns,
         page_margin=c(0, 0.05, 0, 0.05),
         title=my_title, legend=gpar(fontsize=4, col="blue", ncols=1))
  #dev.off()
  
}

#Vectcorlor in function group chosen
colorNbUpDown<-c("plum","plum4",
                 "mediumpurple1","mediumpurple4",
                 "violetred1","violetred4",
                 "LightCyan","Aquamarine","skyblue2","royalblue3",
                 "LightPink","indianred1","firebrick1","red3" )
#color for sankeyNetwork
colorNbUpDown2<-c("#D8BFD8","#DA70D6",
                  "#9370DB","#663399",
                  "#FF69B4","#C71585",
                  "#E0FFFF","#7FFFD4","#6495ED","#00008B",
                  "#FFB6C1","#F08080","#DC143C","#8B0000" )
names(colorNbUpDown)<-c("1_UP.1_DOWN" ,"2_UP.2_DOWN",
                        "1_UP.2_DOWN", "1_UP.3_DOWN" ,
                        "2_UP.1_DOWN","3_UP.1_DOWN",
                        "1_DOWN","2_DOWN","3_DOWN" ,"4_DOWN",
                        "1_UP.","2_UP.","3_UP.","4_UP.")
colorNbconnection<-c("white","#FABDA6","#D3DDDC","#C7B19C","#446455")
names(colorNbconnection)<-c("0","1","2","3","4")

colorCombiDayDEG<-c("#fef6be","#e2f4df","#dbeaf5","#d8d9ea",
                    "#d6efa6","#c2d9ed","#ffcce4","#a0dbb5","#d5b7d7","#8895c2",
                    "#008349","#065fa3","#d42857","#7e2175","#67081d")
names(colorCombiDayDEG)<-c("D0","D2","D4","D7"
                           ,"D0.D2","D0.D4","D0.D7","D2.D4","D2.D7","D4.D7",
                           "D0.D2.D4","D0.D2.D7","D0.D4.D7","D2.D4.D7","D0.D2.D4.D7")


#test for 

tableSankeyplotTOT<-tableSankeyplot(signiDEgene,Typecellv,1,F)

orderTypecell=c("ECs","FAPs","MuSCs","Inflammatory-Mac","Resolving-Mac")
colorsType=c("#10b387ff","#3d85c6ff","#b171f1ff","#ff9900ff","#cc0000ff")
names(colorsType)=orderTypecell
for(typecell in orderTypecell){

pdf(paste0(analyse,"test_All_",typecell,"_sankeyColorUPDOWN.pdf"), width=10, height=7) # Set up PDF device
  
  

tableSankeyplotFAPs<-tableSankeyplot(signiDEgene,typecell,1,T)
sankeyPantaRhei(tableSankeyplotFAPs,colorsType[[typecell]],"NbUpDown",colorNbUpDown,paste0(typecell,"_combisupp1"),F,T)

sankeyPantaRhei(tableSankeyplotFAPs,colorsType[[typecell]],"NbUpDown",colorNbUpDown,paste0(typecell,"_combisupp1"),T,T)
sankeyPantaRhei(tableSankeyplotFAPs,colorsType[[typecell]],"Nbconnection",colorNbconnection,paste0(typecell,"_combisupp1"),F,T)
sankeyPantaRhei(tableSankeyplotFAPs,colorsType[[typecell]],"Nbconnection",colorNbconnection,paste0(typecell,"_combisupp1"),T,T)
sankeyPantaRhei(tableSankeyplotFAPs,colorsType[[typecell]],"CombiDayDEG",colorCombiDayDEG,paste0(typecell,"_combisupp1"),F,T)
sankeyPantaRhei(tableSankeyplotFAPs,colorsType[[typecell]],"CombiDayDEG",colorCombiDayDEG,paste0(typecell,"_combisupp1"),T,T)

tableSankeyplotFAPs<-tableSankeyplot(signiDEgene,typecell,1,F)
sankeyPantaRhei(tableSankeyplotFAPs,colorsType[[typecell]],"NbUpDown",colorNbUpDown,paste0(typecell,"_combisupp1"),F,F)
sankeyPantaRhei(tableSankeyplotFAPs,colorsType[[typecell]],"NbUpDown",colorNbUpDown,paste0(typecell,"_combisupp1"),T,F)
sankeyPantaRhei(tableSankeyplotFAPs,colorsType[[typecell]],"Nbconnection",colorNbconnection,paste0(typecell,"_combisupp1"),F,F)
sankeyPantaRhei(tableSankeyplotFAPs,colorsType[[typecell]],"Nbconnection",colorNbconnection,paste0(typecell,"_combisupp1"),T,F)
sankeyPantaRhei(tableSankeyplotFAPs,colorsType[[typecell]],"CombiDayDEG",colorCombiDayDEG,paste0(typecell,"_combisupp1"),F,F)
sankeyPantaRhei(tableSankeyplotFAPs,colorsType[[typecell]],"CombiDayDEG",colorCombiDayDEG,paste0(typecell,"_combisupp1"),T,F)

tableSankeyplotFAPs<-tableSankeyplot(signiDEgene,typecell,2,T)
sankeyPantaRhei(tableSankeyplotFAPs,colorsType[[typecell]],"NbUpDown",colorNbUpDown,paste0(typecell,"_combisupp2"),F,T)
sankeyPantaRhei(tableSankeyplotFAPs,colorsType[[typecell]],"NbUpDown",colorNbUpDown,paste0(typecell,"_combisupp2"),T,T)
sankeyPantaRhei(tableSankeyplotFAPs,colorsType[[typecell]],"Nbconnection",colorNbconnection,paste0(typecell,"_combisupp2"),F,T)
sankeyPantaRhei(tableSankeyplotFAPs,colorsType[[typecell]],"Nbconnection",colorNbconnection,paste0(typecell,"_combisupp2"),T,T)
sankeyPantaRhei(tableSankeyplotFAPs,colorsType[[typecell]],"CombiDayDEG",colorCombiDayDEG,paste0(typecell,"_combisupp2"),F,T)
sankeyPantaRhei(tableSankeyplotFAPs,colorsType[[typecell]],"CombiDayDEG",colorCombiDayDEG,paste0(typecell,"_combisupp2"),T,T)

tableSankeyplotFAPs<-tableSankeyplot(signiDEgene,typecell,2,F)
sankeyPantaRhei(tableSankeyplotFAPs,colorsType[[typecell]],"NbUpDown",colorNbUpDown,paste0(typecell,"_combisupp2"),F,F)
sankeyPantaRhei(tableSankeyplotFAPs,colorsType[[typecell]],"NbUpDown",colorNbUpDown,paste0(typecell,"_combisupp2"),T,F)
sankeyPantaRhei(tableSankeyplotFAPs,colorsType[[typecell]],"Nbconnection",colorNbconnection,paste0(typecell,"_combisupp2"),F,F)
sankeyPantaRhei(tableSankeyplotFAPs,colorsType[[typecell]],"Nbconnection",colorNbconnection,paste0(typecell,"_combisupp2"),T,F)
sankeyPantaRhei(tableSankeyplotFAPs,colorsType[[typecell]],"CombiDayDEG",colorCombiDayDEG,paste0(typecell,"_combisupp2"),F,F)
sankeyPantaRhei(tableSankeyplotFAPs,colorsType[[typecell]],"CombiDayDEG",colorCombiDayDEG,paste0(typecell,"_combisupp2"),T,F)


dev.off()
}
