#  Explains why GSEA was NOT performed in Tau selected genes for this contrast (Old vs Young by celltype by day)
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
consensusfile <- "Tau/conseTau_ensemblid.rds"
taudefile <- "TauPlusDEinfo_bd_bct.rds"
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
      labs(title = d) +  theme(legend.position = 'none', 
                              axis.text.y=element_text(size=6),
                              axis.text.x=element_text(size=7))
  }
  pdf(paste0(odir, "GSEA_bd_bct/","PROBLEM_tau_DE_cross_bd_bct.pdf"), width=8, height=8)
  title = ggdraw() + draw_label("by day + by celltype,  DE+Tau info")
  captionall = ggdraw() + draw_label(
    paste("PROBLEM detected:\n", 
          "few available Tau for DEGs (absLFC >1) for this contrast\n ",
          "we see that with this contrast, we cannot use Tau, because most of\n",
          "the genes in full DE output  list cannot be associated to a Tau index\n ,", 
          "and the few that could be associated to their Tau exhibit no effect !!!!"),
    size=5
  )
  print( plot_grid(title, plot_grid(plotlist = p_ , ncol=2), legend, 
                   captionall, nrow=4, 
            rel_heights = c(1,10, 1 , 2)) ) # print for pdf inside function OBLIGED
  dev.off()
}  
#makebarplots(mixedfull)
# ======================= END barplots Tau available or not ===================

# we see that with this contrast, we cannot use Tau, because most of
# the genes in full DE output  list cannot be associated to a Tau index , 
# and the few that could be associated to their Tau exhibit no effect !!!!

# ============================  scatter plots conclusion ===================

maketauavailnotsuit <- function(mixedfull, CUTOFFTAU){
  pdf(paste0(odir,"GSEA_bd_bct/", "tau_DE_notSuitableGSEA.pdf"), width = 8, height = 8)
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


# ============================================================================
# approach making UNIQUE symbols by Tau indexes
# ============================================================================
# extract only Tau > CUTOFFTAU and deduplicate with max tau
dedup <- list()
for (d in daysv){
  tf <- consensustau[[d]] %>% filter(Tau >= CUTOFFTAU) %>% 
    group_by(id) %>% slice_max(Tau)  ### ??? is this ok ????
  dedup[[d]] <- tf
} # end for

cat("read how 'consensus' is determined and why, in folder 'exam_INTER_conditions/static/readmedocs' :",
    "\n","file recap_explain.md, section **Tau classification and consensus** ")
# !!  use here function importdedup
 
# ================ concatenate Tau and *full* DEG information ==================
# This uses DEG that resulted from script exam_Inter_cond_sta.R
doconcat = TRUE
if (doconcat){
  # cross information with DEG results FULL version
  degfull <- readRDS(paste0(odir, "rds/shot_rds_full.rds"))
  daysv = c("D0","D2","D4","D7")
  mdegtau <- list()
  for (d in daysv){
    Ddeg <- degfull %>% filter(day == d) %>% select(symbol,log2FoldChange,padj,type)
    
    df_tmp <- right_join(dedup[[d]] %>% mutate(type = whichMAX),
                         Ddeg %>% filter(padj <= 0.5), # this also excludes NA padj
                         by=c("symbol","type"))
    
    df_tmp <- df_tmp %>% filter(abs(log2FoldChange) > 0.1)
    df_tmp <- df_tmp %>% mutate(tauAvailable=ifelse(is.na(Tau), "no", "yes"))
    df_tmp <- df_tmp %>% mutate(Tau_discrete = case_when(
      # is.na(Tau) ~ "NA", # was not available
      Tau >= CUTOFFTAU ~ "SPECIFIC",
      TRUE ~ "HOUSEKEEPING"
    ), LFC_discrete = case_when(
      log2FoldChange >= 1 ~ "A. UP",
      log2FoldChange < 1 & log2FoldChange >= .5  ~ "B. 0.5 to 1",
      log2FoldChange <= -.5 & log2FoldChange > -1 ~ "D. -1 to -0.5",
      log2FoldChange <= -1 ~  "E. DOWN",
      TRUE ~ "C. -0.5 < lfc < 0.5"     )  ) # end mutate
    
    df_tmp <- df_tmp %>% mutate(padj_discrete =  case_when(
      is.na(padj) ~ "NApadj",
      padj <= 0.005 ~ "a. SIGNIFICANT (<= 0.005)",
      padj <= 0.05 & padj > 0.005 ~ "b. SIGNIFICANT (0.05 to 0.005)",
      padj > 0.05 & padj <= 0.5 ~ "c. 0.05 to 0.5",
      padj > 0.5 ~ "d. more than 0.5") )
    df_tmp$day <- d
    mdegtau[[d]] <- df_tmp
  }
  giant = bind_rows(mdegtau)
  # check EACH DAY: there are still some genes classified as "specific" that appear 
  # *duplicated*  among two cell types
  for (d in daysv){
    f <- giant %>% filter(day==d & Tau_discrete=="SPECIFIC") 
    n_occ <- data.frame(table(f$symbol))
    print("checking each day , printing duplicates")
    print(unique(f$type))
    print(n_occ[n_occ$Freq > 1,])
  }
  # make unique (a symbol cannot be repeated across celltypes in one day )
  # criterion : smallest padj
  mdegtau <- list() # clear to reuse this var
  for (d in daysv){
    f <- giant %>% filter(day==d) %>% group_by(symbol) %>%
      slice_min(padj)
    mdegtau[[d]] <- f
  }
  giant <- bind_rows(mdegtau)
  saveRDS(giant, paste0(odir, "rds/", taudefile))
}

# ===================         prepare for MAplot        =======================

giant <- readRDS(paste0(odir, "rds/", taudefile))
for (d in daysv){
  f <- giant %>% filter(day==d & Tau_discrete=="SPECIFIC") 
  n_occ <- data.frame(table(f$symbol))
  print("checking each day , printing duplicates")
  print(unique(f$type))
  if(dim(n_occ[n_occ$Freq > 1,])[1] == 0){ print("ok for plots (no duplics)")
  }else{print("attention, duplicated symbols across celltypes this day")}
}
# =====================================  end load ==============================

# ==================================== MAplot ===============================
giant <- readRDS(paste0(odir, "rds/", taudefile))
#pdf(paste0(odir, "GSEA_bd_bct/", "TauPlusDEinfo_plot_uniqueid.pdf"), width=10, height=6)
ggplot(giant, aes(x=Tau, y=log2FoldChange, color= Tau >= .5) ) + 
  geom_point(size=.5) + geom_hline(yintercept = c(-1.2,1.2),
                                   linetype = "dashed") +
  scale_color_brewer(palette = "Dark2") + facet_grid(padj<=0.05 ~ day) + 
  labs(title="DE information and Tau index", 
       caption="y axis plot labels: TRUE correspond to padj <=0.05") 
#dev.off()
# ====================================  end MAplot ===========================


# =================================  use Mosaic Plot :=======================
# reorder factors : 
giant$Tau_discrete <- ordered(giant$Tau_discrete, levels = c("SPECIFIC", "HOUSEKEEPING"))
revlevelpj <- rev(levels(as.factor(giant$padj_discrete) ))
giant$padj_discrete <- ordered(giant$padj_discrete, levels = revlevelpj)
#revlevelLFC <- rev(levels(as.factor(giant$LFC_discrete) ))
#giant$LFC_discrete <- ordered(giant$LFC_discrete, levels = revlevelLFC)
modx_labs = sapply(levels(as.factor(giant$padj_discrete)), function(x)
  unlist(str_split(x,"\\."))[1])

pdf(paste0(odir, "GSEA_bd_bct/", "TauPlusDEinfo_mosaic.pdf"), width=12, height=10)
ggplot(giant)  +
  geom_mosaic(aes(x=product(LFC_discrete), 
                  conds = product(Tau_discrete), 
                  fill = padj_discrete )) + 
  scale_fill_manual(values=c( "#404788FF", "#238A8DFF",
                              "#55C667FF"), 
                    name="padj") +
  labs(x= "padj", title="Tau, Log2FC and padj values") + 
  scale_x_productlist("padj_categorical", labels=c("c","b","a", "","","")) +
  theme_mosaic() + theme(axis.text.x = element_text(angle = 90),
                         plot.margin = unit(c(.1,.1,.1,.1), "cm"),
                         plot.title= element_text(hjust=0.5)) + coord_flip() +  
  facet_grid(day~Tau_discrete)
dev.off()

# manual colors viridis --> scale_fill_manual(values=c("#440145FF", "#404788FF", "#238A8DFF",
#                              "#55C667FF", "#FDE725FF")),   
# =================================  another Plot :=======================
# proportion of housekeeping and prop specific along celltypes and days

worko <- giant %>% filter(abs(log2FoldChange) >= 0.5) # select only interesting changes

maxpj = round(max(worko$padj),1)
minlfc = round(min(abs(worko$log2FoldChange)),1)

pdf(paste0(odir, "GSEA_bd_bct/", "toplisted_seeTauValues.pdf"), width=6, height=5)
ggplot(worko, aes(x=type, fill= Tau >= CUTOFFTAU)) + geom_bar() +
   labs(title=paste("Number of top listed genes (padj < ",maxpj,", abs(lfc) >",minlfc,")"),
        caption= paste("NA values mean that no Tau value was found.",
                       "CUTOFFTAU is :",CUTOFFTAU))
dev.off()


# end
