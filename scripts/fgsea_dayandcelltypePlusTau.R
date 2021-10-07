#  Explains why GSEA was not performed in Tau selected genes.
#  for this contrast (Old vs Young by celltype by day)
# --
# johaGL
library(fgsea)
library(tidyverse)
library(msigdbr)
library(cowplot)
setwd("~/BulkAnalysis_plusNetwork/")

CUTOFFTAU = 0.3
daysv = c('D0', 'D2', 'D4', 'D7')
DEfile = "rds/shot_rds_full.rds"
odir = "exam_INTER_conditions/static/"
infoinputfile = "GSEA/csv/infoinput-V2"  # csv or md_txt extensions 
pdfplotfile = "preGSEA_byday_bycelltype-V2.pdf"
gseards = "GSEA/rds/"
consensusfile <- "Tau/conseTau_ensemblid.rds"

consensustau <- readRDS(consensusfile)
DEdf <- readRDS(paste0(odir, DEfile))
summary(DEdf)
# names(DEdf)
# [1] "baseMean"       "log2FoldChange" "lfcSE"          "pvalue"         "padj"           "id"             "day"           
# [8] "type"           "symbol"        
# > names(consensustau)
# [1] "D0" "D2" "D4" "D7"

if (!"symbol" %in% colnames(DEdf)){
  genes_df <- read.table("data/genesinfo.csv",sep="\t",header=T)
  DEdf$symbol <- genes_df[match(DEdf$id, genes_df$Geneid),]$symbol
  head(DEdf)
  print("added symbols to dataframe")
}else{ print("symbol already in dataframe")}

concatenateDEwithTau <- function(DEdf, consensustau){
  mixedfull <- list()
  for (d in names(consensustau)){
    print(paste(" %%%%%%%%%%%%%", d, "%%%%%%%%%%%%%"))
    dem.here <- DEdf %>% filter(day == d) %>% filter(!is.na(padj))
    tau.here <- consensustau[[d]]
    tmp_ <- list()
    for (c in unique(dem.here$type)){
      print(paste("     ", c))
      mix <- left_join(
          dem.here %>% filter(type==c),
          tau.here %>% filter(whichMAX==c),
          by = "id"  )# end _join
      tmp_[[c]] <- mix %>% mutate(tauAvailable=ifelse(is.na(Tau), "no", "yes"))
    } # end for
    mixedfull[[d]] <- bind_rows(tmp_)
  } # end for
  return(mixedfull)
}
  
mixedfull <- concatenateDEwithTau(DEdf, consensustau)
# tail(mixedfull$D2 %>% filter(!is.na(Tau) & !is.na(padj)) %>% arrange(desc(padj)))
# baseMean log2FoldChange lfcSE   pvalue     padj id                 day   type   symbol.x symbol.y Tau               class        whichMAX nbMAX exclusiveOld
# <dbl>          <dbl> <dbl>    <dbl>    <dbl> <chr>              <chr> <chr>  <chr>    <chr>    <chr>             <chr>        <chr>    <chr> <chr>       
#   1     429.           2.02 0.221 6.49e-21 1.35e-17 ENSMUSG00000026768 D2    sCs    Itga8    Itga8    0.767054165861923 intermediate sCs      1     0  

# ============================ barplots Tau available or not ===================
makebarplots <- function(mixedfull){
  p_ <- list()
  legend <- cowplot::get_legend(ggplot(mixedfull$D0, aes(x=type, fill=padj<0.05))+
                  geom_bar() + theme(legend.direction = 'horizontal')) 
  for (d in names(mixedfull)){
    p_[[d]] <- ggplot(data = mixedfull[[d]] %>%
        mutate(tauAvailable=ifelse(is.na(Tau), "TauUnavailable","TauAvailable")), 
           aes(x=type, fill=padj<0.05)) + geom_bar() + 
      facet_grid(abs(log2FoldChange)>= 1 ~ tauAvailable, 
            labeller=labeller(.rows = c('TRUE'="absFC>1", 'FALSE'="absFC<1"))) +
      labs(title = d) + theme(legend.position = 'none', 
                              axis.text.y=element_text(size=6),
                              axis.text.x=element_text(size=7))
  }
  pdf(paste0(odir, "tau_DE_cross_bd_bct.pdf"), width=8, height=8)
  title = ggdraw() + draw_label("by day + by celltype,  DE+Tau info")
  print( plot_grid(title, plot_grid(plotlist = p_ , ncol=2), legend, nrow=3, 
            rel_heights = c(1,10, 1)) ) # print for pdf inside function OBLIGED
  dev.off()
}  
#makebarplots(mixedfull)
# ======================= END barplots Tau available or not ===================

# we see that with this contrast, we cannot use Tau, because most of
# the genes in full DE output  list cannot be associated to a Tau index , 
# and the few that could be associated to their Tau exhibit no effect !!!!

# ============================  scatter plots conclusion ===================

maketauavailnotsuit <- function(mixedfull, CUTOFFTAU){
  pdf(paste0(odir,"tau_DE_notSuitableGSEA.pdf"), width = 8, height = 8)
  agg <- ggplot(bind_rows(mixedfull) %>% mutate(Tau = as.numeric(Tau)) %>% 
                  filter(tauAvailable == "yes"),
                aes(x=Tau, y=as.numeric(log2FoldChange), color=padj<=0.05 )) + 
    geom_point(size=.4, alpha=0.6) + scale_color_manual(values=c("darkgray", "darkgreen")) +
    geom_vline(xintercept = CUTOFFTAU, linetype="dashed", size=.2) +
    facet_grid(type~day) + theme_light() + 
    labs(title="Contrast OvsY by celltype by day, Tau available only", 
         caption = paste("this figure shows that a selection based on Tau",
                         "for this contrast, is NOT suitable for GSEA \n",
                         "specific genes (Tau>0.3) for ECs and macrophages do not change with age"))
  print(agg)
  dev.off()
}
#maketauavailnotsuit(mixedfull, CUTOFFTAU)
# end
