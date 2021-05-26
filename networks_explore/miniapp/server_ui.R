library(shiny)
library(ggplot2)
library(dplyr)
library(igraph)
library(visNetwork)

mywdir <- "~/BulkAnalysis_plusNetwork/networks_explore/"
grdir <- "graphobjs/"
g <- read_graph(paste0(mywdir,grdir,"myminigraph.ml"), format="graphml")
vertex_attr(g)$numid <- vertex_attr(g)$id
edge_attr(g)$color <- edge_attr(g)$ecolor


data <- toVisNetworkData(g)
data$nodes$label = data$nodes$genesym
data$nodes$value = log(data$nodes$averagexp)*10
data$nodes$groupname = data$nodes$celltype
data$edges$width = data$edges$weight*10
#nodes = as.data.frame(vertex_attr(g))
#edges = as.data.frame(edge_attr(g))

u <- shinyUI(fluidPage(
  titlePanel("very small test, by joha 2021"),
  
  sidebarLayout(position = "left",
                sidebarPanel( h4("Parameters"), width = 3,
                              sliderInput("range", "Range:",
                                          min = 1, max = 1000,
                                          value = c(200,500)),
                              sliderInput("JOKER","this joker:",0,10,2,step=0.1),
                              numericInput("ngroups","Groups:",5,1,11,1),
                              actionButton("GO","Generate"),
                              actionButton("SAVE", "SAVE"),
                              textOutput("networkstat")
                ),
                mainPanel(h4("Network Plots"), width = 9,
                          tabsetPanel(
                            tabPanel("visNetwork",
                                     visNetworkOutput("visnetwork",
                                     height="500px"),
                                     style = "background-color: #ffffff;"),
                            tabPanel("otherTag", 
                                     style= "background-color: #eeeeee;" )
                          ),
                          
                )
  )
))

s <- function(input, output){
  
  net <- reactiveValues(nodes=NULL,edges=NULL,groups=NULL)
  
  observeEvent(
   
    input$GO,{
    print("regenerating network")
    ng <- input$ngroups
    JOKER <- input$JOKER
    
    net$nodes <- data$nodes
    net$edges <- data$edges
    
    output$visnetwork <- renderVisNetwork({ 
      req(net$edges)
      rangee <-  input$range
      print(rangee)
      print(unique(net$nodes$groupname))
      netout <- visNetwork(net$nodes,net$edges,
                           physics = FALSE) %>%
        visEdges(shadow=T, arrows="to", physics=FALSE) %>%
        visInteraction(navigationButtons = TRUE) %>%
        visGroups(groupname="M1")
      netout
      
    })# ebd renderVisNetwork
      
    output$networkstat <- renderText({
      sprintf("\nNodes:%d  Edges:%d Groups:%d",
              nrow(net$nodes),nrow(net$edges),nrow(net$groups))
    })
    
  }#end input$GO
  ) # end observeEvent
  
  observe(input$SAVE,
  )
  
  
}

shinyApp(u,s)