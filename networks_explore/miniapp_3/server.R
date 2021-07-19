library(shiny)
library(ggplot2)
library(dplyr)
library(igraph)
library(visNetwork)
library(tidyverse)
library(DT)
library(shinyalert)

#mywdir <- "~/BulkAnalysis_plusNetwork/networks_explore/miniapp_3/"
grdir <- "graphobjs_copy/"
source("ui.R") #TODO:  set good ui
source("serverFun.R")
DEclassic_file <- "DE_copy/shot_dataframe.csv"

# ===================== declare empty reactiveValues ==========================
igElems_list <- reactiveValues(x=list())
currentday <- reactiveValues(x='')
subnet.y <- reactiveValues(nodes=NULL,edges=NULL,groups=NULL)
subnet.o <- reactiveValues(nodes=NULL,edges=NULL,groups=NULL)
DEclassic.show <- reactiveValues(x=NULL)
crossedtab <- reactiveValues(x=NULL)

# ===================== function that returns subgraph ========================
doinduced <- function(g, interestvec, orderinp){
  vertex_attr(g)$numid <- vertex_attr(g)$id
  vertex_attr(g)$uniname <- vertex_attr(g)[["_nx_name"]]
  selnodes = V(g)[uniname %in% interestvec]
  selegoV <- ego(g, order=orderinp, nodes = selnodes, mode = "all", mindist = 0)
  # turn the returned list of igraph.vs objects into a graph
  selegoG <- induced_subgraph(g,unlist(selegoV))
  return(selegoG)
}

# ############################ SERVER function ################################
server <- function(input, output, session ){
  # ====================== Condition to start  ========================
  #print(length(isolate(igElems_list$x)))
  if (length(isolate(igElems_list$x)) == 0 ){
    output$labmain_Young <- renderText({
      paste("Press a DAY, then button 'Load only'")
    })
  }
  # ============================ Event select day  ===========================
  observeEvent(
    input$DAY,{
      print(paste(isolate(currentday$x), "<===chosenday"))
      if((!isolate(currentday$x) == "") & (!isolate(currentday$x)==input$DAY)){
        output$visnet_y <- renderPlot({NULL})
        output$visnet_o <- renderPlot({NULL})
        output$labmain_Young <- renderText({"....please load..."})
        output$labmain_Old <- renderText({".."})
      }
    } # end input$DAY
  ) # end observeEvent
  # ===================== Event load tables (renderDT)   =======================
  observeEvent(
    input$LOADONLY,{  # button load
      print("you pressed button to load the igs")
      currentday$x <- isolate(input$DAY)
      for (age in c("Young","Old")){
        tmpReacGraph = reactive({
          inFile = currentday$x
          if (!is.null(inFile))
            read.graph(file=paste0(grdir, age,"_",currentday$x,
                                  "_igraph_unfi.ml"), format="graphml")
        })
        igElems_list$x[[age]] <- isolate(tmpReacGraph())
        print(paste("okloaded", age, currentday$x))
      } # end for age in c("Young","Old")
      # ---------------------- inside Event load tables,....  ------------------------
      output$labmain_Young <- renderText({
       paste("Young, ", currentday$x, "loaded internally, \n go to tables to see available nodes")})
      output$labmain_Old <- renderText({
       paste("Old, ", currentday$x, "loaded internally, \n go to tables to see available nodes")})
      #return(currentday)
      output$visnet_y <- renderPlot({NULL})
      output$visnet_o <- renderPlot({NULL})
      ###  yield tables
      output$tableyoung <- renderDT(
        igraph::as_data_frame(igElems_list$x[["Young"]],"vertices"),
        filter = "top",
        options = list(pageLength=15))
        # renderTable({ igraph::as_data_frame(igElems_list$x[["Young"]],"vertices")})
      output$tableold <- renderDT(
        igraph::as_data_frame(igElems_list$x[["Old"]], "vertices"),
        filter = "top",
        options = list(pageLength=15)) 
     } # end input LOADONLY
   ) # end observeEvent LOADONLY
     
  # ================== Event show entire graphs (heavy)   ====================
  observeEvent(
    input$SHOWMAIN,{
      if (currentday == isolate(currentday$x)){
        print("ok to show")
      }else{ output$labmain_Young <- renderText({"YOU HAVEN'T LOADED ANYTHING !!!! "})}
      print(currentday)
      print(paste(isolate(currentday$x),"====> is the current day selected"))
      print(isolate(input$DAY))
      output$visnet_y <- renderPlot({
        coords = layout_with_kk(igElems_list$x[["Young"]])
        plot.igraph(igElems_list$x[["Young"]],vertex.size=5, mark.border=NA,
                    edge.border=NA,  vertex.border="white",
                    vertex.label =NA, edge.arrow.size=0.005,
                    edge.arrow.width=0.00001,
                    layout=coords)
      })
      output$visnet_o <- renderPlot({
        coords = layout_with_kk(igElems_list$x[["Old"]])
        plot.igraph(igElems_list$x[["Old"]],vertex.size=5, mark.border=NA,
                    edge.border=NA,  vertex.border="white",
                    vertex.label =NA, edge.arrow.size=0.005,
                    edge.arrow.width=0.00001,
                    layout=coords)
      })
      # ------------------- put titles
      output$labmain_Young <- renderText({
        paste("Young, ", currentday$x)
      })
      output$labmain_Old <- renderText({
        paste("Old, ", currentday$x)
      })
      # ------------------- send genes list ... pending
      print("ok show all (very slow)")
    }) # end observeEvent SHOWMAIN
  # ====================== Crossing old young tables  =======================
  observeEvent(
    input$docross,{
      print(length(igElems_list$x))
      if (length(igElems_list$x) == 0){
        output$crossinf <- renderText({"error: no networks were loaded !!"})
        crossedtab$x = NULL
      }else{
        V_vec = list()
        for (age in c("Young","Old")){
          V_vec[[age]] = igraph::as_data_frame(igElems_list$x[[age]], "vertices")[["_nx_name"]]
          print(head(igraph::as_data_frame(igElems_list$x[[age]], "vertices"), 2))
        }
        crossedtab$x[["ExclusiveYoung"]] = sort(setdiff(V_vec[["Young"]],V_vec [["Old"]]))
        crossedtab$x[["ExclusiveOld"]] = sort(setdiff(V_vec[["Old"]],V_vec [["Young"]]))
        crossedtab$x[["Intersection"]] = sort(intersect(V_vec[["Old"]], V_vec[["Young"]])) 
        # output$crossinf <- renderText({"chevere"})
        output$crossinf <- renderPrint(crossedtab$x)
        
      } #end else
    } #end input$docross
  )
  
  # ====================== Event showing subgraphs  ========================
  observeEvent(
    input$GO,{
      massageDATA <- function(myg){
        data <- toVisNetworkData(myg)
        data$nodes$label = data$nodes$genesym
        data$nodes$value = data$nodes$specificity * 10  # !! $averagexp changed!
        data$nodes$groupname = data$nodes$celltype
        data$edges$width = data$edges$weight * 10
        return(data)
      }
      if (currentday$x == ''){
        shinyalert("Oops!", "You did not loaded a DAY", type = "error")
        print(paste("!!!", currentday$x))
      }else{
        pickneigh = list("Young"=c(''), "Old"=c(''))
        here.subgr = list()
        # ---------------------- young visnet ---------------------------------
        if (input$Young_nodes == "" ){
          output$young <- renderVisNetwork({NULL})
          output$labtextyoung <- renderText({
            paste("NO nodes declared to search for. Introduce a single one:\n",
                  "Il34_FAPs\n or a comma separated list:\n ",
                  "Crlf1_sCs, Edil3_sCs")})
        }else{
          pickneigh[["Young"]] = unlist( str_split(
            gsub(" ","",input$Young_nodes), ",") )
          ytrycatch <- tryCatch({
            here.subgr[["Young"]] <- doinduced(
              igElems_list$x[["Young"]], pickneigh[["Young"]], orderinp = input$NEIGH )
            print(head(igraph::as_data_frame(here.subgr[["Young"]], "vertices"),1))
            subnet.y$nodes <- massageDATA(here.subgr[["Young"]])$nodes
            subnet.y$edges <- massageDATA(here.subgr[["Young"]])$edges
            print("ok subnet young")
            output$labtextyoung <- renderText({paste("YOUNG day ", currentday$x)})
            output$young <- renderVisNetwork({
              req(subnet.y$edges)
              print(unique(subnet.y$nodes$groupname))
              
              netout <- visNetwork(subnet.y$nodes,
                                   subnet.y$edges,
                                   physics = FALSE) %>%
                visEdges(shadow=T, arrows="to", physics=FALSE) %>%
                visInteraction(navigationButtons = TRUE)
              netout
            }) # end renderVisnetwork
          }, warning = function(war){
            output$young <- renderVisNetwork({NULL})
            output$labtextyoung <- renderText({
              paste("this node is not not in Young at this day")})
            print(war)
          }, error = function(err){
            output$young <- renderVisNetwork({NULL})
            output$labtextyoung <- renderText({
              paste("this node is not not in Young at this day")})
            print(err)
          }, finally = { print("bad entry by user at input$Young_nodes")}
          ) # end ytrycatch
        } # end else young
        # ---------------------- old visnet ---------------------------------
        if (input$Old_nodes == "" ){
          output$old = renderVisNetwork({NULL})
          output$labtextold <- renderText({
            paste("NO nodes declared to search for. Introduce a single one:\n",
                  "Lgi3_sCs\n or a comma separated list:\n",
                  "Wnt3_FAPs, Kremen1_ECs")})
        }else{
          pickneigh[["Old"]] = unlist( str_split(
            gsub(" ", "",input$Old_nodes), ",") )
          
          otrycatch <- tryCatch({
            here.subgr[["Old"]] <- doinduced(
              igElems_list$x[["Old"]],
              pickneigh[["Old"]],
              orderinp = input$NEIGH )
            subnet.o$nodes <- massageDATA(here.subgr[["Old"]])$nodes
            subnet.o$edges <- massageDATA(here.subgr[["Old"]])$edges
            print("ok subnet old")
            output$labtextold <- renderText({paste("OLD, day ", currentday$x)})
            output$old <- renderVisNetwork({
              req(subnet.o$edges)
              print(unique(subnet.o$nodes$groupname))
              netout <- visNetwork(subnet.o$nodes,
                                   subnet.o$edges,
                                   physics = FALSE) %>%
                visEdges(shadow=T, arrows="to", physics=FALSE) %>%
                visInteraction(navigationButtons = TRUE)
              netout
            }) #end renderVisNetwork
            }, warning = function(war) {
              output$old = renderVisNetwork({NULL})
              output$labtextold <- renderText({
                paste("no this node is not in Old at this day")})
            }, error = function(err){
              output$old = renderVisNetwork({NULL})
              output$labtextold <- renderText({
                paste("no this node is not in Old at this day")})
            }, finally = { print("bad param by user in input$Old_nodes ")}
            )#end otrycatch
          } # end else old
        } # end if else for shinyalert when currentday$x == ''
      }#end input$GO
    ) # end observeEvent GO  (button Animated)
  # =========================================================================
  # end showing subgraphs
  # ========================== Event showing DE =============================
  observeEvent(
    input$DEclassical,{
      DEGs <- read.table(DEclassic_file, sep=',',header=T) 
      DEclassic.show$x <- DEGs %>% select(day,type, symbol, id, log2FoldChange, padj, baseMean) 
      rm(DEGs)
      output$DE_classical <- renderDT(
        DEclassic.show$x,
        filter = "top",
        options = list(pageLength=50, widthChange= TRUE,
                       columnDefs = list(list(visible=FALSE, targets=c('baseMean'))))
      )#end renderDT
    } # end input$DEclassical
  )
  
  # ====================== Event showing intersection DE LR ===================
  observeEvent(
    input$execu_LR_DE,{
        if (currentday$x == ''){
          output$daywaspick <- renderText({paste("you did not loaded DAY in L-R section")})
          print(isolate(input$DAY))
        }else if(currentday$x != isolate(input$DAY)){
          output$daywaspick <- renderText({paste("error, the day in selector is not same as loaded")})
        }else{
          output$daywaspick <- renderText({paste("the DAY you picked is : ", currentday$x)})
          # subset DE by current day and set 'Symbol_celltype' in "name" column
          tmp_c <- DEclassic.show$x %>% filter(day == currentday$x)
          tmp_c <- tmp_c %>% mutate(name = paste0(symbol, '_',type))
          print(head(tmp_c),2)
          v_c_youngUp <- tmp_c %>% filter(log2FoldChange > 0) %>% pull(name) 
          v_c_oldUp <- tmp_c %>% filter(log2FoldChange < 0) %>% pull(name) 
          print(v_c_youngUp)
          print(v_c_oldUp)
          # get both networks old and young and do intersections
          g.o = igraph::as_data_frame(igElems_list$x[["Old"]], "vertices")
          g.y = igraph::as_data_frame(igElems_list$x[["Young"]], "vertices")
          v_o = g.o[["_nx_name"]]
          v_y = g.y[["_nx_name"]]
          union_net = union(v_o, v_y)
          #print(inDEG_notinNET)
          intersectDEG_n_NET = intersect(union_net, tmp_c$name)
          print(paste("** ==========", currentday$x,"=================== **"))
          print("** checking for picked day DEG genes and comparing with Network **")
          #print(head(v_c_youngUp))
          #print(head(v_y))
          YinDEG_notinNET = setdiff(v_c_youngUp, v_y)
          OinDEG_notinNET = setdiff(v_c_oldUp, v_o)
          inNET_notinDEG_hyperYoung = setdiff(crossedtab$x[["ExclusiveYoung"]], v_c_youngUp )
          inNET_notinDEG_hyperOld = setdiff(crossedtab$x[["ExclusiveOld"]], v_c_oldUp)
          print(paste(length(intersectDEG_n_NET), "<== intersection DEG and NET"))
          print(paste(length(YinDEG_notinNET), "<== hyperYoungDEG not in NET"))
          print(paste(length(OinDEG_notinNET), "<== hyperOldDEG not in NET"))
          print(paste(length(inNET_notinDEG_hyperYoung), "<== exclusiveYoung but not seen in DEG list"))
          print(paste(length(inNET_notinDEG_hyperOld), "<== exclusiveOld but not seen in DEG list"))
          print("======================================================================")
          #print(paste(length(, " ????"))
        }
      }
    )# end event LR_DE
  } # end server function



