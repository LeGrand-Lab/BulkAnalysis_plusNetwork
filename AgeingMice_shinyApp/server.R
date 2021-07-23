# johaGL 2021
library(shiny)
library(ggplot2)
library(dplyr)
library(igraph)
library(visNetwork)
library(tidyverse)
library(DT)
library(shinyalert)
library(circlize)
library(kableExtra)
library(heatmap3)

# TODO: !!somewhere take note: 'dynamic' DE, for M1 yielded no significant 
#  because from one day to the next the change was not important
# that is why I take only old vs youn in staticsnapshot comparisons

mywdir <- "~/BulkAnalysis_plusNetwork/AgeingMice_shinyApp/"
setwd(mywdir)
grdir <- "graphobjs_copy/"
source("ui.R") #TODO:  set good ui
DEclassic_file <- "DE_copy/shot_dataframe.csv"
aggreg_matrices <- readRDS("Data/aggreg_matrices.rds")
mhp = readRDS("Data/matrices4heatmaps.rds")

# ===================== declare empty reactiveValues ==========================
currentday <- reactiveValues(x='')
igElems_list <- reactiveValues(x=list())  # the networks objects, both ages
aggreg_day <- reactiveValues(x=list())
subnet.y <- reactiveValues(nodes=NULL,edges=NULL,groups=NULL)
subnet.o <- reactiveValues(nodes=NULL,edges=NULL,groups=NULL)
DEclassic.show <- reactiveValues(x=NULL)
crossedtab <- reactiveValues(x=NULL)

# global variable colors
cellcolors = list(
  "ECs"="#0072B2",
  "FAPs"="#F0E442",
  "M1" = "#D55E00",
  "M2" =  "#CC79A7",
  "Neutro" =  "#009E73",
  "sCs" = "#56B4E9" )

# ================== function that returns chorddiagram =======================
dochord <- function(matrx, unicolors){
  chordDiagram(matrx, directional = 1,grid.col = unlist(cellcolors),
               col = unicolors,
               direction.type = c("diffHeight", "arrows"),
               link.arr.type = "big.arrow",
             self.link = 1,
             annotationTrack = c("name", "grid"," axis"))} 

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
  reacAggregPlots <- reactiveValues(doPlot = FALSE)

  # ====================== Condition to start  ========================
  #print(length(isolate(igElems_list$x)))
  if (length(isolate(igElems_list$x)) == 0 ){
    output$labmain_young <- renderText({
      paste("Select a DAY, then button 'Load only'")
    })
  }
  # ========================== display DEGs by default =========================
  displayDEGS <- function(){
    DEGs <- read.table(DEclassic_file, sep=',',header=T) 
    DEclassic.show$x <- DEGs %>% select(day,type, symbol, id, log2FoldChange, padj, baseMean) 
    rm(DEGs)
    output$DE_classical <- renderDT(
      DEclassic.show$x,
      filter = "top",
      options = list(pageLength=50, widthChange= TRUE,
                     columnDefs = list(list(visible=FALSE, targets=c('baseMean'))))
    )
  }
  displayDEGS()
  # ============================ Event select day  ===========================
  observeEvent(
    input$DAY,{
      print(paste(isolate(currentday$x), "<===chosenday"))
      if((isolate(currentday$x) == "") || (isolate(currentday$x)!= input$DAY)){
        reacAggregPlots$doPlot <- FALSE
        output$radiores_young_title <- renderText({""})
        output$radiores_young_A <- renderPlot({NULL})
        output$radiores_young_B <- renderTable({NULL})
        output$radiores_old_title <- renderText({""})
        output$radiores_old_A <- renderPlot({NULL})
        output$radiores_old_B <- renderTable({NULL})
        output$labmain_young <- renderText({"....please load..."})
        output$labmain_old <- renderText({".."})
        output$crossinf <- renderText({NULL})
      }else{
       
        next
      }
    } # end input$DAY
  ) # end observeEvent
  
  # ===================== Event load tables (intro plots (chord) and renderDT  =======================

  # REQUIRES  function 'displayintroBEST'   :: strong test
  displayintroBEST <- function(ag){
    print(names(ag))  # Young part:
    output$radiores_young_title <- renderText({
      if (reacAggregPlots$doPlot == FALSE) { return () } else {
        if (input$radioAggreg == "Ratio"){
          return("Cumulated weight / number of connections ")
        }else{
          return(isolate(input$radioAggreg))} # endif
      } }) # endif # end renderText
    output$radiores_young_A <-  renderPlot({
      if (reacAggregPlots$doPlot == FALSE) { return(NULL) } else {
        dochord(ag[["Young"]][[input$radioAggreg]], unicolors= ag[["Colorsmat"]])
      } } , height = 450, width = 450 )
    output$radiores_young_B <-  function(){  # the table is special
      if (reacAggregPlots$doPlot == FALSE) { return () } else {
        df <- ag[["Young"]][[input$radioAggreg]]
        df %>%
          knitr::kable("html") %>%
          kable_styling("striped", full_width = F) %>%
          add_header_above(c(" ", "To" = dim(df)[2])) %>% 
          add_footnote("rownames are 'From' nodes, colnames are 'To' nodes")
      } } #end table # end young part
    # for old : 
    output$radiores_old_title <- renderText({
      if (reacAggregPlots$doPlot == FALSE) { return () } else {
        if (input$radioAggreg == "Ratio"){
          return("Cumulated weight / number of connections ")
        }else{
          return(isolate(input$radioAggreg))
        }
      } })
    output$radiores_old_A <-  renderPlot({
      if (reacAggregPlots$doPlot == FALSE) { return(NULL) } else {
        dochord(ag[["Old"]][[input$radioAggreg]], unicolors= ag[["Colorsmat"]])
      } } , height = 450, width = 450 )
    output$radiores_old_B <-  function(){  # the table is special
      if (reacAggregPlots$doPlot == FALSE) { return () } else {
        df <- ag[["Old"]][[input$radioAggreg]]
        df %>%
          knitr::kable("html") %>%
          kable_styling("striped", full_width = F) %>%
          add_header_above(c(" ", "To" = dim(df)[2])) %>% 
          add_footnote("rownames are 'From' nodes, colnames are 'To' nodes")
      } } #end table # end old part
  } #end function displayintroBEST
  
  observeEvent(
    input$LOADSTART,{  # button load
      # by default 
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
        aggreg_day$x[[age]] <- aggreg_matrices[[age]][[currentday$x]]
        print("ok aggreg day from aggreg_matrices")
      } # end for age in c("Young","Old")
      reacAggregPlots$doPlot <- TRUE
      print(reacAggregPlots$doPlot)
      output$thisiscondition <- renderText({reacAggregPlots$doPlot})
     
      if (reacAggregPlots$doPlot == TRUE){
        displayintroBEST(aggreg_day$x)
        output$labmain_young <- renderText({
          paste("Young, ", currentday$x, "loaded, go to tables to see available nodes")})
        
        output$labmain_old <- renderText({
          paste("Old, ", currentday$x, "loaded, go to tables to see available nodes")})
        # ----------------------   load DT tables,....  ------------------------
        ###  yield tables
        output$labtabyoung <- renderText({paste(currentday$x, " ! ")})
        output$tableyoung <- renderDT(
          igraph::as_data_frame(igElems_list$x[["Young"]],"vertices"),
          filter = "top",
          options = list(pageLength=15))
        # renderTable({ igraph::as_data_frame(igElems_list$x[["Young"]],"vertices")})
        output$labtabold <- renderText({paste(currentday$x, " ! ")})
        output$tableold <- renderDT(
          igraph::as_data_frame(igElems_list$x[["Old"]], "vertices"),
          filter = "top",
          options = list(pageLength=15)) 
      }
     } # end input LOADSTART
   ) # end observeEvent LOADSTART
  
  
  # ============================ reactive radioButton displays ===========================
  ## show aggregated visuals (chorddiagrams and 'From-To' matrices)
  observeEvent(
    input$radioAggreg,{
      if (reacAggregPlots$doPlot == FALSE){
        print("YOU HAVNT LOADED Anything")
              }else{
        displayintroBEST(aggreg_day$x)
      }
    }
  )
  #            
 # end radioButtons
 
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
        if (input$young_nodes == "" ){
          output$young <- renderVisNetwork({NULL})
          output$labtextyoung <- renderText({
            paste("NO nodes declared to search for. Introduce a single one:\n",
                  "Il34_FAPs\n or a comma separated list:\n ",
                  "Crlf1_sCs, Edil3_sCs")})
        }else{
          pickneigh[["Young"]] = unlist( str_split(
            gsub(" ","",input$young_nodes), ",") )
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
          }, finally = { print("bad entry by user at input$young_nodes")}
          ) # end ytrycatch
        } # end else young
        # ---------------------- old visnet ---------------------------------
        if (input$old_nodes == "" ){
          output$old = renderVisNetwork({NULL})
          output$labtextold <- renderText({
            paste("NO nodes declared to search for. Introduce a single one:\n",
                  "Lgi3_sCs\n or a comma separated list:\n",
                  "Wnt3_FAPs, Kremen1_ECs")})
        }else{
          pickneigh[["Old"]] = unlist( str_split(
            gsub(" ", "",input$old_nodes), ",") )
          
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
            }, finally = { print("bad param by user in input$old_nodes ")}
            )#end otrycatch
          } # end else old
        } # end if else for shinyalert when currentday$x == ''
      }#end input$GO
    ) # end observeEvent GO  (button Animated)
  # =========================================================================
  # end showing subgraphs
  
  # ====================== Event showing intersection DE LR ===================
  observeEvent(
    input$execu_LR_DE,{
        if (currentday$x == ''){
          output$daywaspick <- renderText({paste("you did not loaded DAY in L-R section")})
          print(isolate(input$DAY))
        }else if(currentday$x != isolate(input$DAY)){
          output$daywaspick <- renderText({paste("error, the day in selector is not same as loaded")})
        }else if(is.null(crossedtab$x)){
          output$daywaspick <- renderText({paste("day:",currentday$x, 
                                                 " you must also click 'docross' in 'cross' tab in L-R section")})
        }else{
          output$daywaspick <- renderText({paste("the DAY you picked is : ", currentday$x)})
          # subset DE by current day and set 'Symbol_celltype' in "name" column
          head(DEclassic.show$x)
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
          intersectDEG_n_NET = intersect(union_net, tmp_c$name)
          YinDEG_notinNET = setdiff(v_c_youngUp, v_y)
          OinDEG_notinNET = setdiff(v_c_oldUp, v_o)
          inNET_notinDEG_hyperYoung = setdiff(crossedtab$x[["ExclusiveYoung"]], v_c_youngUp )
          inNET_notinDEG_hyperOld = setdiff(crossedtab$x[["ExclusiveOld"]], v_c_oldUp)
         
          print(paste("** %%%%%%%%%%%%%", currentday$x,"%%%%%%%%%%%%% **"))
          print("** checking for picked day DEG genes and comparing with Network **")
          print(paste(length(intersectDEG_n_NET), "<== intersection DEG and NET"))
          print(paste(length(YinDEG_notinNET), "<== hyperYoungDEG not in NET"))
          print(paste(length(OinDEG_notinNET), "<== hyperOldDEG not in NET"))
          print(paste(length(inNET_notinDEG_hyperYoung), "<== exclusiveYoung but not seen in DEG list"))
          print(paste(length(inNET_notinDEG_hyperOld), "<== exclusiveOld but not seen in DEG list"))
          print("%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
          
          output$details_crossing <- renderTable({
            data.frame("description" = c("intersection DEG and NET",
                                         "UpregulatedYoungDEG not in NET",
                                         "UpreglatedOldDEG not in NET"
                                         ),
                       "genesymbols_celltype" = c(
                         paste(intersectDEG_n_NET, collapse=", "),
                         paste(YinDEG_notinNET, collapse=", " ),
                         paste(OinDEG_notinNET, collapse=", ")
                         
                       ),
                       "number" = c(length(intersectDEG_n_NET),
                                    length(OinDEG_notinNET ),
                                    length(OinDEG_notinNET)
                       )
            ) # end data.frame
            })
        }
      }
    )# end event LR_DE
  
  
  # ========================== addd pathways (display by default) ============================
  sapply(names(mhp[['D0']]), function(CT, k = 'D0'){
    print(CT)
    output[[paste0(k,"_",CT,"_","UP")]] <- renderPlot({
      heatmap3(as.matrix(mhp[[k]][[CT]][["UP"]]),  
               scale = "none",
               col=colorRampPalette(c("gray","firebrick3"))(256), 
               Colv = NA, Rowv = NA,
               cexRow = 1,
               cexCol = 1,
               ColSideWidth = ncol(mhp[[k]][[CT]][["UP"]]),
               margins = c(7,10),
               main = paste(k, CT, "UP"))
    })
    output[[paste0(k,"_",CT,"_","DOWN")]] <- renderPlot({
      heatmap3(as.matrix(mhp[[k]][[CT]][["DOWN"]]),  
               scale = "none",
               col=colorRampPalette(c("navy","gray"))(256), 
               Colv = NA, Rowv = NA,
               cexRow = 1,
               cexCol = 1,
               ColSideWidth = ncol(mhp[[k]][[CT]][["DOWN"]]),
               margins = c(7,10),
               main = paste(k, CT, "DOWN"))
    })
  }) # end sapply D0
  sapply(names(mhp[['D2']]), function(CT, k = "D2"){
    print(CT)
    output[[paste0(k,"_",CT,"_","UP")]] <- renderPlot({
      heatmap3(as.matrix(mhp[[k]][[CT]][["UP"]]),  
               scale = "none",
               col=colorRampPalette(c("gray","firebrick3"))(256), 
               Colv = NA, Rowv = NA,
               cexRow = 1,
               cexCol = 1,
               ColSideWidth = ncol(mhp[[k]][[CT]][["UP"]]),
               margins = c(7,10),
               main = paste(k, CT, "UP"))
    })
    output[[paste0(k,"_",CT,"_","DOWN")]] <- renderPlot({
      heatmap3(as.matrix(mhp[[k]][[CT]][["DOWN"]]),  
               scale = "none",
               col=colorRampPalette(c("navy","gray"))(256), 
               Colv = NA, Rowv = NA,
               cexRow = 1,
               cexCol = 1,
               ColSideWidth = ncol(mhp[[k]][[CT]][["DOWN"]]),
               margins = c(7,10),
               main = paste(k, CT, "DOWN"))
    })
  }) # end sapply D2
  sapply(names(mhp[['D4']]), function(CT, k='D4'){
    print(CT)
    output[[paste0(k,"_",CT,"_","UP")]] <- renderPlot({
      heatmap3(as.matrix(mhp[[k]][[CT]][["UP"]]),  
               scale = "none",
               col=colorRampPalette(c("gray","firebrick3"))(256), 
               Colv = NA, Rowv = NA,
               cexRow = 1,
               cexCol = 1,
               ColSideWidth = ncol(mhp[[k]][[CT]][["UP"]]),
               margins = c(7,10),
               main = paste(k, CT, "UP"))
    })
    output[[paste0(k,"_",CT,"_","DOWN")]] <- renderPlot({
      heatmap3(as.matrix(mhp[[k]][[CT]][["DOWN"]]),  
               scale = "none",
               col=colorRampPalette(c("navy","gray"))(256), 
               Colv = NA, Rowv = NA,
               cexRow = 1,
               cexCol = 1,
               ColSideWidth = ncol(mhp[[k]][[CT]][["DOWN"]]),
               margins = c(7,10),
               main = paste(k, CT, "DOWN"))
    })
  }) # end sapply D4
  sapply(names(mhp[['D7']]), function(CT, k='D7'){
    print(CT)
    output[[paste0(k,"_",CT,"_","UP")]] <- renderPlot({
      heatmap3(as.matrix(mhp[[k]][[CT]][["UP"]]),  
               scale = "none",
               col=colorRampPalette(c("gray","firebrick3"))(256), 
               Colv = NA, Rowv = NA,
               cexRow = 1,
               cexCol = 1,
               ColSideWidth = ncol(mhp[[k]][[CT]][["UP"]]),
               margins = c(7,10),
               main = paste(k, CT, "UP"))
    })
    output[[paste0(k,"_",CT,"_","DOWN")]] <- renderPlot({
      heatmap3(as.matrix(mhp[[k]][[CT]][["DOWN"]]),  
               scale = "none",
               col=colorRampPalette(c("navy","gray"))(256), 
               Colv = NA, Rowv = NA,
               cexRow = 1,
               cexCol = 1,
               ColSideWidth = ncol(mhp[[k]][[CT]][["DOWN"]]),
               margins = c(7,10),
               main = paste(k, CT, "DOWN"))
    })
  }) # end sapply D7
  # ========================== end added pathways ==============================
  } # end server function




#### END
# NOTE :  # universal colors used for cell types
# ok for color blindness!  hex codes: https://rdrr.io/cran/ggthemes/man/colorblind.html
