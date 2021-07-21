# ===================
# just parts to add if necessary

# !! hey:  "output$labmain_..." and "output$visnet_..." must be changed by new stuff
# not yet created
# intact event as before:
# ================== Event show entire graphs (heavy)   ====================

observeEvent(
  input$SHOWMAIN,{
    if (currentday == isolate(currentday$x)){  # BAD CONDITION !!! CHECK
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

################################### done this morning 21 july:
# ============================ reactive radioButton displays
## show aggregated visuals (chorddiagrams and 'From-To' matrices)
observeEvent(
  input$radioAggreg,{
    df <- head(mtcars)
    if (reacAggregPlots$doPlot == FALSE){
      print("YOU HAVNT LOADED Anything")
    }else{
      print(input$radioAggreg)
      output$radiores_young_title <- renderText({input$radioAggreg})
      output$radiores_young_A <- renderPlot({plot(c(1,0), c(0,1))})
      output$radiores_young_B <- function(){
        df %>%
          knitr::kable("html") %>%
          kable_styling("striped", full_width = F) %>%
          add_header_above(c(" ", "To" = dim(df)[2])) %>% 
          add_footnote("rownames are 'From' nodes, colnames are 'To' nodes")
      }
    }
  }
)
#             kable_styling("striped", full_width = F) %>%




