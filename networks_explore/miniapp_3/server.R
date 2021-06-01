library(shiny)
library(ggplot2)
library(dplyr)
library(igraph)
library(visNetwork)
library(tidyverse)
mywdir <- "~/BulkAnalysis_plusNetwork/networks_explore/"
grdir <- "graphobjs/"
source(paste0(mywdir,"miniapp_3/ui.R"))

myigobs <- reactiveValues(x=list())
currentday <- reactiveValues(x='')
subnet.o <- reactiveValues(nodes=NULL,edges=NULL,groups=NULL)
subnet.y <- reactiveValues(nodes=NULL,edges=NULL,groups=NULL)

mywdir <- "~/BulkAnalysis_plusNetwork/networks_explore/"
grdir <- "graphobjs/"

doinduced <- function(g, interestvec, orderinp){
  vertex_attr(g)$numid <- vertex_attr(g)$id
  vertex_attr(g)$uniname <- vertex_attr(g)[["_nx_name"]]
  selnodes = V(g)[uniname %in% interestvec]
  selegoV <- ego(g, order=orderinp, nodes = selnodes, mode = "all", mindist = 0)
  # turn the returned list of igraph.vs objects into a graph
  selegoG <- induced_subgraph(g,unlist(selegoV))
  return(selegoG)
}
# 
# g <- read_graph(paste0(mywdir,grdir,"Young_D2_igraph.ml"), format="graphml")
# #g <- read_graph(paste0(mywdir,grdir,"myminigraph.ml"), format="graphml")
# vertex_attr(g)$numid <- vertex_attr(g)$id
# edge_attr(g)$color <- edge_attr(g)$ecolor
# 
# data <- toVisNetworkData(g)
# data$nodes$label = data$nodes$genesym
# data$nodes$value = log2(data$nodes$averagexp)*10
# data$nodes$groupname = data$nodes$celltype
# data$edges$width = data$edges$weight*10
# #nodes = as.data.frame(vertex_attr(g))
# #edges = as.data.frame(edge_attr(g))
# MAXWEIGHT.EDGE = max(data$edges$width)

server <- function(input, output){
  #print(length(isolate(myigobs$x)))
  if (length(isolate(myigobs$x)) == 0 ){
    output$init_young <- renderText({
      paste("Press a DAY, then button Load ")
    })
  }
  observeEvent(
    input$DAY,{
      print(isolate(currentday$x))
      if((!isolate(currentday$x) == "") & (!isolate(currentday$x)==input$DAY)){
        output$testigraph <- renderPlot({NULL})
        output$testigraph_b <- renderPlot({NULL})
        output$init_young <- renderText({"....please load..."})
        output$init_old <- renderText({".."})
      }
    }
  )
  
   # =================== button Load 
   observeEvent(
     input$LOADONLY,{
       currentday$x <- isolate(input$DAY)
       for (age in c("Young","Old")){
         tmpReacGraph = reactive({
           inFile = currentday$x
           if (!is.null(inFile))
             read.graph(file=paste0(mywdir,grdir, age,"_",currentday$x,
                                    "_igraph.ml"), format="graphml")
         })
         myigobs$x[[age]] <- isolate(tmpReacGraph())
         print(paste("okloaded", age, currentday))
       } # end for age in c("Young","Old")
       
       output$init_young <- renderText({
         paste("Young, ", currentday$x, "loaded internally, press Show")})
       output$init_old <- renderText({
         paste("Old, ", currentday$x, "loaded internally, press Show")})
       #return(currentday)
       output$testigraph <- renderPlot({NULL})
       output$testigraph_b <- renderPlot({NULL})
     }
   ) # end observeEvent LOADONLY
     
  # ===================== button Show 
  observeEvent(
    input$SHOWMAIN,{
      print(currentday)
      print(paste(isolate(currentday$x),"hihihih"))
      print(isolate(input$DAY))
      print("haha")
      
      output$testigraph <- renderPlot({
        coords = layout_with_kk(myigobs$x[["Young"]])
        plot.igraph(myigobs$x[["Young"]],vertex.size=5, mark.border=NA,
                    edge.border=NA,  vertex.border="white",
                    vertex.label =NA, edge.arrow.size=0.005,
                    edge.arrow.width=0.00001,
                    layout=coords)
      })
      output$testigraph_b <- renderPlot({
        coords = layout_with_kk(myigobs$x[["Old"]])
        plot.igraph(myigobs$x[["Old"]],vertex.size=5, mark.border=NA,
                    edge.border=NA,  vertex.border="white",
                    vertex.label =NA, edge.arrow.size=0.005,
                    edge.arrow.width=0.00001,
                    layout=coords)
      })
      # ------------------- put titles
      output$init_young <- renderText({
        paste("Young, ", currentday$x)
      })
      output$init_old <- renderText({
        paste("Old, ", currentday$x)
      })
      # ------------------- send genes list ... pending
      print("ok show all (very slow)")
    }) # end observeEvent SHOWMAIN
  
  # =================== button Animated
  observeEvent(
     input$GO,{
       if (input$Young_nodes == "" & input$Old_nodes==""){
         output$labtextyoung <- renderText({
           "MUST ! declare nodes to select"})
       }else{
         pickneigh = list( "Young"= unlist( str_split(
                            gsub(" ","",input$Young_nodes),
                            ",") ),
                          "Old"= unlist( str_split(
                            gsub(" ", "",input$Old_nodes),
                            ",") ) ) # end list
         here.subgr = list()
         # =================== get induced subgraph
         for(age in c("Young","Old")){
          here.subgr[[age]] <- doinduced(
            myigobs$x[[age]],
            pickneigh[[age]],
            orderinp = input$NEIGH
          )
         }
         print("ok til here 1 ")
         massageDATA <- function(myg){
           data <- toVisNetworkData(myg)
           data$nodes$label = data$nodes$genesym
           data$nodes$value = data$nodes$averagexp * 10
           data$nodes$groupname = data$nodes$celltype
           data$edges$width = data$edges$weight * 10
           return(data)
         }
         print(head(igraph::as_data_frame(here.subgr[["Young"]], "vertices"),1))
         print(head(igraph::as_data_frame(here.subgr[["Old"]], "vertices"),1))
         print("ok , an eye on both induced sg ")
         # --------------------- young visnet
         if (dim(igraph::as_data_frame(here.subgr[["Young"]], "vertices"))[1] == 0){
           print(paste("no graph to return in YOUNG : ", pickneigh[["Young"]]))
         }else{
           subnet.y$nodes <- massageDATA(here.subgr[["Young"]])$nodes
           subnet.y$edges <- massageDATA(here.subgr[["Young"]])$edges
           print("ok subnet young")
           
           output$labtextyoung <- renderText({"YOUNG day ...."})
           output$young <- renderVisNetwork({
             req(subnet.y$edges)
             print(unique(subnet.y$nodes$groupname))
             netout <- visNetwork(subnet.y$nodes,
                                  subnet.y$edges,
                                  physics = FALSE) %>%
               visEdges(shadow=T, arrows="to", physics=FALSE) %>%
               visInteraction(navigationButtons = TRUE) 
             netout
           })# ebd renderVisNetwork young
         } # end ifelse young
         # -------------------- old visnet  
         if (dim(igraph::as_data_frame(here.subgr[["Old"]], "vertices"))[1] == 0){
           print(paste("no graph to return in OLD : ",pickneigh[["Old"]]))
         }else{
           subnet.o$nodes <- massageDATA(here.subgr[["Old"]])$nodes
           subnet.o$edges <- massageDATA(here.subgr[["Old"]])$edges
           print("ok subnet old")
           
           output$labtextold <- renderText({"OLD, day ...."})
           output$old <- renderVisNetwork({
             req(subnet.o$edges)
             print(unique(subnet.o$nodes$groupname))
             netout <- visNetwork(subnet.o$nodes,
                                  subnet.o$edges,
                                  physics = FALSE) %>%
               visEdges(shadow=T, arrows="to", physics=FALSE) %>%
               visInteraction(navigationButtons = TRUE) 
             netout
           })# ebd renderVisNetwork old
         } # end ifelse old
       } # end else if input$..._nodes is filled
      }#end input$GO
   ) # end observeEvent GO  (button Animated)
  } # end server function

#shinyApp(ui,server)


# #  select induced graph
# #View(igraph::as_data_frame(g, "edges"))
# #View(igraph::as_data_frame(g,"vertices"))
# interest = c("Dhh_sCs")
# selegoG = doinduced(g, interest, 7)
# # plot the subgraph
# plot(selegoG, vertex.size=5, mark.border=NA,
#      edge.border=NA,  vertex.border="white",
#      vertex.label =NA, edge.arrow.size=0.005,
#      edge.arrow.width=0.00001)
# # https://ourednik.info/maps/2018/09/21/create-a-subgraph-from-the-neighborhood-of-specific-vertices-in-igraph/
# 
# 
# # ----------------------- 
# #print(typeof(myigobs$x[["Young"]]))  #  !!! animated !!!
# #output$labtextyoung <- renderText({input$Young_nodes})
# #output$labtextold <- renderText({input$Old_nodes})


# 
#   # # for the moment only loading "old" data
#   # net$nodes <- data$nodes
#   # net$edges <- data$edges
#   # 
#   # # FAKKE !!! young data:
#   # net_y$nodes <- data$nodes %>% filter(averagexp < 30)
#   # net_y$edges <- data$edges %>% filter(from %in% net_y$nodes$id)
#   # print(net_y$edges)
#   # 
#   # print("!!!!!")
#   # print(input$range)
#   # 
#   # output$labtextyoung <- renderText({"YOUNG"})
#   # output$young <- renderVisNetwork({ 
#   #   req(net$edges)
#   #   rangee <-  input$range
#   #   print(rangee)
#   #   print(unique(net$nodes$groupname))
#   #   netout <- visNetwork(net$nodes,net$edges,
#   #                        physics = FALSE) %>%
#   #     visEdges(shadow=T, arrows="to", physics=FALSE) %>%
#   #     visInteraction(navigationButtons = TRUE) %>%
#   #     visGroups(groupname="M2")
#   #   netout
#   #   
#   # })# ebd renderVisNetwork 
#   # 
#   # output$labtextold <- renderText({"OLD"})
#   # output$old <- renderVisNetwork({ 
#   #   req(net_y$edges)
#   #   rangee <-  input$range
#   #   print(rangee)
#   #   netout <- visNetwork(net_y$nodes,net_y$edges,
#   #                        physics = FALSE) %>%
#   #     visEdges(shadow=T, arrows="to", physics=FALSE) %>%
#   #     visInteraction(navigationButtons = TRUE) %>%
#   #     visGroups(groupname="M1")
#   #   netout
#   #   
#   # })# ebd renderVisNetwork old
#   
#   #output$networkstat <- renderText({
#    # sprintf("\nNodes:%d  Edges:%d Groups:%d",
#     #        nrow(net$nodes),nrow(net$edges),nrow(net$groups))
#   #})
#   
