library(shiny)
library(ggplot2)
library(dplyr)
library(igraph)
library(visNetwork)

locgr <- "~/BulkAnalysis_plusNetwork/graphmodeling/"
g <- read_graph(paste0(locgr,"myminigraph.ml"), format="graphml")
vertex_attr(g)$numid <- vertex_attr(g)$id
edge_attr(g)$color <- edge_attr(g)$ecolor


data <- toVisNetworkData(g)
data$nodes$label = data$nodes$genesym
data$nodes$value = log2(data$nodes$averagexp)*10
data$nodes$groupname = data$nodes$celltype
data$edges$width = data$edges$weight*10
#nodes = as.data.frame(vertex_attr(g))
#edges = as.data.frame(edge_attr(g))
MAXWEIGHT.EDGE = max(data$edges$width)

# tryig to make shadows on groups by vertex cell type
v <- as_data_frame(g, "vertices")
v$numid.num<- sapply(v$id, function(x) stringr::str_replace(x,"n", ""))
plot(g,
     mark.groups = list(
       unlist(v %>% filter(celltype=="M2" ) %>% select(numid.num)),
       unlist(v %>% filter(celltype=="M1" )  %>% select(numid.num)),
       unlist(v %>% filter(celltype=="FAPs")  %>% select(numid.num)),
       unlist(v %>% filter(celltype=="sCs") %>% select(numid.num))
     ),mark.col=c("#C5E5E7","#ECD89A","#C5E5E7",
                  "#ECD89A","#E66101"), mark.border=NA)
####  Error in simple_vs_index(x, ii, na_ok) : Unknown vertex selected


ui <- shinyUI(fluidPage(
  titlePanel("very small test, by joha 2021"),
  
  sidebarLayout(position = "left",
                sidebarPanel( h4("Parameters"), width = 3,
                              sliderInput("range", "Range:",
                                          min = 1, max = 1000,
                                          value = c(200,500)),
                              sliderInput("JOKER","this joker:",0,10,2,step=0.1),
                              numericInput("ngroups","Groups:",5,1,11,1),
                              actionButton("SAVE", "SAVE"),
                              actionButton("GO","Generate animated"),
                              textOutput("networkstat")
                ),
                mainPanel(h4("Network Plots"), width = 9,
                          tabsetPanel(
                            tabPanel("Main",
                                     #visNetworkOutput("Main", ="350px"),
                                     plotOutput("testigraph"),
                                     plotOutput("testigraph_b"),
                                     style = "background-color: #ffffff;"),
                            tabPanel("animated", 
                                     verbatimTextOutput("haha"),
                                     visNetworkOutput("young",
                                                      height="350px"),
                                     verbatimTextOutput("hoho"),
                                     visNetworkOutput("old",
                                                      height="350px"),
                                     style= "background-color: #ffffff;" )
                          ),
                          
                )
  )
))

server <- function(input, output){
  
  output$testigraph <- renderPlot({
    
    gg <- make_ring(10)
    values <- lapply(1:10, function(x) sample(1:10,3))
    if (interactive()) {
      plot(gg, vertex.shape="pie", vertex.pie=values,
           vertex.pie.color=list(heat.colors(5)),
           vertex.size=seq(10,30,length=10), vertex.label=NA)
    } # end if interactive
    } )
  output$testigraph_b <- renderPlot({
    
    plot.igraph(g, layout=layout.mds)
  })
  net <- reactiveValues(nodes=NULL,edges=NULL,groups=NULL)
  net_y <- reactiveValues(nodes=NULL,edges=NULL,groups=NULL)
  
  observeEvent(
   
    input$GO,{
    print("regenerating network")
    ng <- input$ngroups
    JOKER <- input$JOKER
    
    # for the moment only loading "old" data
    net$nodes <- data$nodes
    net$edges <- data$edges
    
    # FAKKE !!! young data:
    net_y$nodes <- data$nodes %>% filter(averagexp < 30)
    net_y$edges <- data$edges %>% filter(from %in% net_y$nodes$id)
    print(net_y$edges)
    
    print("!!!!!")
    print(input$range)
    
    output$haha <- renderText({"OLD"})
    output$young <- renderVisNetwork({ 
      req(net$edges)
      rangee <-  input$range
      print(rangee)
      print(unique(net$nodes$groupname))
      netout <- visNetwork(net$nodes,net$edges,
                           physics = FALSE) %>%
        visEdges(shadow=T, arrows="to", physics=FALSE) %>%
        visInteraction(navigationButtons = TRUE) %>%
        visGroups(groupname="M2")
      netout
      
    })# ebd renderVisNetwork 
    
    output$hoho <- renderText({"YOUNG"})
    output$old <- renderVisNetwork({ 
      req(net_y$edges)
      rangee <-  input$range
      print(rangee)
      netout <- visNetwork(net_y$nodes,net_y$edges,
                           physics = FALSE) %>%
        visEdges(shadow=T, arrows="to", physics=FALSE) %>%
        visInteraction(navigationButtons = TRUE) %>%
        visGroups(groupname="M1")
      netout
      
    })# ebd renderVisNetwork old
    
    output$networkstat <- renderText({
      sprintf("\nNodes:%d  Edges:%d Groups:%d",
              nrow(net$nodes),nrow(net$edges),nrow(net$groups))
    })
    
  }#end input$GO
  ) # end observeEvent
}

shinyApp(ui,server)