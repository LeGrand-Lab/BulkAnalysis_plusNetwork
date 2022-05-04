# Change manually to get ALL, day and type plot


GSEAtempo<-GSEAsigni %>% mutate( pathway2 = str_replace_all(pathway, "REACTOME_", "") )
prefix_name="ALL"
t=""
GSEAtempo<-GSEAsigni %>% mutate( pathway2 = str_replace_all(pathway, "REACTOME_", "") )
  
listpathway<- c(unique(GSEAtempo$pathway2),"NOT_IN_REACTOME_HIERARCHIE","IN_REACTOME_HIERARCHIE")
  
data_nested_2<-data.frame()
sizeVector<-c()
dico_path<-data.frame()
dico_path<-rbind(dico_path,c("IN_REACTOME_HIERARCHIE",rep(0,13)))
colnames(dico_path)<-c("path",1:13)
  pathwayFound<-c()
  
for(l in rownames(data_nested)){
    newrow=rep("", 13)
    sizetempo=0
    pathFound=F
    for (c in 13:1){
      if (is.na(data_nested[l,c]!= "") == F){
        if ( data_nested[l,c] %in% listpathway){
          if (pathFound == F){
            sizetempo=GSEAtempo %>% filter(pathway2 == data_nested_full[l,c] ) %>% select( size ) %>% unlist() %>% as.character() %>% as.numeric()
            sizetempo=sum(sizetempo,na.rm = T )
          }
          pathFound=T 
          newrow[[c]]<- str_replace_all(data_nested[l,c],'_',' ')
        }
        else{
          if (pathFound == T){
            if (data_nested[l,c] %in% dico_path$path){
              newrow[[c]]<- paste0(c,'_',dico_path[dico_path$path == data_nested[l,c], c+1])
            }
            else {
              newKey<-dico_path[length(dico_path$path),]
              newKey[["path"]]<-data_nested[l,c]
              newKey[[c+1]]<-as.numeric(newKey[[c+1]])+1
              newrow[[c]]<- paste0(c,'_',newKey[[c+1]])
              dico_path<-rbind(dico_path,newKey)
            }
          }
        }
      }
    }
    if (pathFound == T){
      newrowExist<-t(data.frame(newrow))
      colnames(newrowExist)<-colnames(data_nested_2)
      if (nrow(merge(newrowExist,data_nested_2))==0){
        data_nested_2<-rbind(data_nested_2,newrow)
        sizeVector<-c(sizeVector,unique(sizetempo))
      }
    }
  }
  
  colnames(data_nested_2)<-c("level1","level2","level3","level4","level5","level6","level7","level8","level9","level10","level11","level12","level13")
  data_nested_2$value=sizeVector
  
  if (dir.exists(paste0(odir,HierarchieData,prefix_name,"/")) == F) {
    dir.create(paste0(odir,HierarchieData,prefix_name,"/"))
  }
  write.table(data_nested_2, file=paste0( odir,HierarchieData,prefix_name,"/",prefix_name,"_table_REACTOME_Hierachie.txt"))
  
  
  data_nested_2$pathString <- paste("roots", data_nested_2$level1, data_nested_2$level2, data_nested_2$level3, data_nested_2$level4,data_nested_2$level5,data_nested_2$level6,data_nested_2$level7, data_nested_2$level8, data_nested_2$level9, data_nested_2$level10,data_nested_2$level11,data_nested_2$level12, sep = "/")
  
  
  data_Node <-as.Node(data_nested_2)
  if (t ==""){
    p <- circlepackeR(data_Node, size = "value",color_min = paste0("hsl(56,80%,80%)"), color_max = paste0("hsl(341,30%,30%)"), width = 1800, height = 1000)
    
  }else{
    p <- circlepackeR(data_Node, size = "value",color_min = paste0("hsl(",hsl_color[[t]],",80%,80%)"), color_max = paste0("hsl(",hsl_color[[t]],",30%,30%)"), width = 1500, height = 1000)
  }
  p
  
  saveWidget(p, file=paste0( odir,plotpathwaysHierarchie,prefix_name,"_REACTOME_Hierachie.html"))


