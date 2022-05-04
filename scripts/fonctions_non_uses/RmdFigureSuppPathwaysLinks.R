      if (length(dataFrameTotalUPDOWN %>% filter(day==d,type==t) %>% subset(symbol %in% genes) %>% select(log2FoldChange) %>% unlist() %>% as.numeric()) != length(genes)){
        tempo=dataFrameTotalUPDOWN %>% filter(day==d,type==t) %>% subset(symbol %in% genes)
        fullGSEA.tempo <- data.frame(tempo ,genes_df[match(tempo$id, genes_df$Geneid),]$symbol)
        extractGeneDuplicate <- fullGSEA.tempo[duplicated(tempo$symbol),]$symbol
        genes2 <-genes[ !genes %in%  extractGeneDuplicate]
        
        genes<-data_frame("symbol"=genes2) %>% as.data.frame()
        
        matrix<-data.frame(symbol=genes$symbol )
        dataFrameTotalUPDOWN %>% filter(day==d,type==t) %>% subset(symbol %in%  genes$symbol) %>% select(symbol) %>% unlist()
        
        matrix <- matrix %>% mutate("random" = dataFrameTotalUPDOWN %>% filter(day==d,type==t) %>% subset(symbol %in%  genes$symbol) %>% select(log2FoldChange) %>% unlist() %>% as.numeric())
      }
      else{
      matrix <- matrix %>% mutate("random" = dataFrameTotalUPDOWN %>% filter(day==d,type==t) %>% subset(symbol %in% genes$symbol) %>% select(log2FoldChange) %>% unlist() %>% as.numeric() )
      }
      colnames(matrix)<-colNAmesMatrix
      
      
      ## Dynamic of log2foldchange genes DE in this pathway through days and  cell type for young versus old 

--------------

```{r DynamiueGenesPathways, echo=FALSE,warning=FALSE,message=FALSE, fig.show='asis'}

PlotGenelogFacrossDay2 <- ggplot(params$tablePlottest, aes(x=day, y=log2F, group=paste(symbol,type,age))) +
      scale_y_continuous(limits=c(params$minlogF,params$maxlogF),) + # fill=name allow to automatically dedicate a color for each group
      geom_line() +
      scale_color_gradient2(midpoint=0,  low="blue", mid="grey88",high="red", limits = c(params$minlogF,params$maxlogF),name="log2FoldChange")+
      geom_hline(yintercept = 0, color="black", linetype="dashed",lwd=0.5,slope=F)+
      geom_point(aes(col=log2F,size=log10(meanCountNormalized+1)))+
      scale_size(name="Counts normalised",limits=c(1,params$maxCountNormalized),breaks = c(log10(10),log10(100),log10(1000),log10(50000)),labels=c('10',"100","1000","50000"))+
      new_scale("color") +
      facet_grid(vars(type),vars(age),
                  labeller=labeller(type=params$TypeUnique))+
      geom_text_repel(data=subset2,
                      aes(label=symbol,color=type),  size=3,min.segment.length =0, segment.size = .8,point.padding= 0.8, force=5, force_pull = 0.1, max.overlaps=15,show.legend=F)  +
      scale_color_manual(values = params$col2, name="Top 10 genes padj<0.05")+
      ylab("log2FoldChange")+ggtitle(unique(params$tablePlottest$Pathway))+
      theme_linedraw()+ theme(strip.text  = element_text(size = 12 , face = "bold"),strip.background = element_rect(fill="grey"),panel.grid.major.x = element_line(color = "gray60", size = 0.8),axis.title= element_text(face="bold"),axis.text = element_text(face="bold",size=14))
    print(PlotGenelogFacrossDay2)
```

