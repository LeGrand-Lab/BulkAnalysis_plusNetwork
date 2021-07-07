library(shiny)
library(ggplot2)
library(dplyr)
library(igraph)
library(visNetwork)
library(DT)

ui <- shinyUI(fluidPage(
  titlePanel("Skeletal Muscle: Ligand Receptor Networks"),
  sidebarLayout(position = "left",
                sidebarPanel( h4("Parameters"), width = 3,
                              selectInput("DAY","DAY ( = timepoint):", list('D0','D2','D4','D7'),
                                          selected = NULL),
                              actionButton("LOADONLY", "Load only (internal)"),
                              actionButton("SHOWMAIN", "Show full nets (slow)"),
                              textInput("Young_nodes","Young network, nodes to select", value=""),
                              textInput("Old_nodes","Old network, nodes to select", value=""),
                              numericInput("NEIGH","neighbors desired", value=1,
                                           min = 1,
                                           max = 50),
                              actionButton("GO","Generate animated"),
                              sliderInput("rangeEdges", "Range edges specificities (weight):",
                                          min = 0.01, max = 1,
                                          value = c(0.1,1)),
                              sliderInput("rangeNodes", "Range vertices (nodes) specificities:",
                                          min = 0.01, max = 1,
                                          value = c(0.1,1)),
                              #sliderInput("JOKER","a pourvoir:",0,10,2,step=0.1),
                              actionButton("SAVE", "SAVE"),
                              textOutput("networkstat")
                ),
                mainPanel(h4("Network Plots"), width = 9,
                          tabsetPanel(
                            tabPanel("Main",
                                     #visNetworkOutput("Main", ="350px"),
                                     verbatimTextOutput("labmain_Young"),
                                     plotOutput("testigraph"),
                                     verbatimTextOutput("labmain_Old"),
                                     plotOutput("testigraph_b"),
                                     style = "background-color: #ffffff;"),
                            tabPanel("table Young", 
                                     verbatimTextOutput("labtabyoung"),
                                     #tableOutput("tableyoung"),
                                     DTOutput("tableyoung"),
                                     style= "background-color: #ffffff;" ),
                            tabPanel("table Old", 
                                     verbatimTextOutput("labtabold"),
                                     DTOutput("tableold"),
                                     style= "background-color: #ffffff;" ),
                            tabPanel("animated", 
                                     verbatimTextOutput("labtextyoung"),
                                     visNetworkOutput("young",
                                                      height="350px"),
                                     verbatimTextOutput("labtextold"),
                                     visNetworkOutput("old",
                                                      height="350px"),
                                     style= "background-color: #ffffff;" )
                       
                          
                ), # end tabsetpanel
                          
        ) # end mainpanel
  ) # end sidebarLayout
))
