library(shiny)
library(ggplot2)
library(dplyr)
library(igraph)
library(visNetwork)
library(DT)
library(shinythemes)
library(shinyjs)

stylebtn1 <- "color: #ffffff; background-color: #288ba8; border-color: #288ba8;"
# stylebtn1 is Mango Tango  .... orange like

# note: nested tabPanels REQUIRE  'sidebarLayout' as outer nest ! 

ui <- navbarPage("Muscle and Age",
# ============================== L-R pairs tabPanel ============================               
  tabPanel( "L-R pairs" ,
     sidebarLayout(
        sidebarPanel( h4("Parameters"), width = 3,
                 selectInput("DAY","DAY (or timepoint):", list('D0','D2','D4','D7'),
                             selected = NULL),
                 actionButton("LOADONLY", "Load only",
                              class="btn-primary"),
                 actionButton("SHOWMAIN", "Full nets (slow)"),
                 br(),
                 br(), 
                 textInput("Young_nodes","Young network, nodes to select", value=""),
                 textInput("Old_nodes","Old network, nodes to select", value=""),
                 numericInput("NEIGH","neighbors desired", value=1,
                              min = 1,
                              max = 50),
                 actionButton("GO","Generate animated", class="btn-success"),
                 br(),
                 sliderInput("rangeEdges", "Range edges specificities (weight):",
                             min = 0.01, max = 1, value = c(0.1,1)),
                 sliderInput("rangeNodes", "Range vertices (nodes) specificities:",
                             min = 0.01, max = 1, value = c(0.1,1)),
                 #sliderInput("JOKER","a pourvoir:",0,10,2,step=0.1),
                 actionButton("SAVE", "SAVE"),
                 textOutput("networkstat")
               ) , # end sidebarPanel

               mainPanel(
                 h2("Ligand-Receptor networks"), width = 9,
                    tabsetPanel(
                         tabPanel("Load",
                                  #visNetworkOutput("Main", ="350px"),
                                  strong(textOutput("labmain_Young")),
                                  plotOutput("testigraph"),
                                  strong(textOutput("labmain_Old")),
                                  plotOutput("testigraph_b"),
                                  style = "background-color: #ffffff;"),
                          
                          tabPanel("table Young", 
                                   textOutput("labtabyoung"),
                                   #tableOutput("tableyoung"),
                                   DTOutput("tableyoung"),
                                   style= "background-color: #ffffff;" ),
                         
                          tabPanel("table Old", 
                                   textOutput("labtabold"),
                                   DTOutput("tableold"),
                                   style= "background-color: #ffffff;" ),
                         
                         tabPanel("cross",
                                  actionButton("docross", "docrossing"),
                                  verbatimTextOutput("crossinf")),
                         
                         tabPanel("animated",
                                  textOutput("labtextyoung"),
                                  visNetworkOutput("young",
                                                   height="350px"),
                                  textOutput("labtextold"),
                                  visNetworkOutput("old",
                                                   height="350px"),
                                  style= "background-color: #ffffff;" )
                    ) # end tabsetpanel
                  ) # end mainPanel
              ) #end sidebarLayout 
           ),
  # ============================== DE tabPanel ============================               

  tabPanel("DE ... in construction"),  
  
  # ============================== About ============================           
  tabPanel("About",
           fluidRow(
            column(2),
             column(5, 
                    h2("About us"),
                    br(),
                    p(paste("\nThis application contains bioinformatic analysis",
                            "related to Bulk RNA-seq experiment conducted on skeletal muscle of",
                            "young and ageing backgrounds in mice. Main authors are",
                            "H. Dieu,  ,  ....., W. Jarassier,  F. Le Grand, B. Chazaud, and  -- Feigé ")),
                    br(),
                    p(paste("Bioinformatic analysis presented here (including this shiny", 
                            "application you are exploring right now)",
                            "was performed by Johanna Galvis, M2 Intern in Bioinformatics 2021")),
                    style="background-color : #f1f0fa"
                    ),
            column(1), #empty column
            column(2, h4("Links to labs"),
                         style="background-color: #b2abd2")
           ) #end fluidRow
        
  )
           
)
