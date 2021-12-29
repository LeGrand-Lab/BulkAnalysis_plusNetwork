##  closely related to 'calc_Tau_Specificity.R'
# this script is only for figure purposes related to Tau index
#  for preprocessing if 'plotsprep_calc_needed' set to TRUE
#  for ulterior visuals, if 'ulterior' set to TRUE
## JohaGL 2021

library(dplyr)
library(tidyverse)
library(openxlsx)
library(reshape2)
library(VennDiagram)
library(ggplot2)
library(ggthemes)
library(cowplot)
require(gridExtra) # calling grid.arrange

setwd("~/BulkAnalysis_plusNetwork/")
plotsprep_calc_needed = F # iff need re-run preprocessing figures
ulterior = T
resdir = 'Tau/'# place where 'calc_Tau_Specificity.R' has saved Tau tables 

# minimal example by hand tau specificity from Yanai:
vec <- c(0, 8, 0, 0, 0, 2, 0, 2, 0, 0, 0, 0)
num <- sum(1-(vec/max(vec)))
den <- length(vec)-1
num/den

calculateTau <- function(vec){
  if(max(vec) == 0){ # avoid zero division error
    tau.index <- 0
  }else{  
    tau.index <- (sum(1-(vec/max(vec))))/(length(vec)-1)
  }
  return(tau.index)
}


# print tau  preprocessing plots if needed,
if (plotsprep_calc_needed){
  age =  "Young"  # age = "Old"
  pdf(paste0("Tau/","plotDetailscalc_", age,".pdf"))
  days=c("D0","D2","D4","D7")
  par(mfrow = c(4,3) )
  for (i in days){
    ktab <- read.table(paste0("data/meanTPM_",age,i,".txt"), sep='\t', header=T,
                       row.names=1) 
    logtab <- log10(ktab+1) 
    hist(unlist(logtab),xlim=c(-3,8), col="wheat", prob=T, cex.main=.9,
         main=paste("compare logarithmic transform\n (TPM+1)",i, age), 
         xlab = "bars:log10 (lines:log2) TPM+1")
    lines(density(unlist(log(ktab+1))),  col="red", cex.main=.7)
    text(4,0.15,"log base 2", col="red")
    # get rid homogeneously low rows in matrix, setting a min quantile cutoff
    myquantile = 0.25
    q.cutoff <- quantile(unlist(logtab),myquantile)
    keep <- apply(logtab, 1, function(x) sum(x >= 0) == length(x) &
                    sum(x > q.cutoff) >= 1)  #and at least one over this cutoff
    logtab <- logtab[keep,] 
    hist(unlist(logtab), col="dodgerblue", 
         main=paste("After filter: retain rows log10(TPM+1)\n",
                    "having at least one value > quantile", myquantile, "(",round(q.cutoff,4),")"  ),
         cex.main=.7, xlab = "log10(TPM+1)")
    print(dim(logtab))
    tau_res <- tibble("id"=rownames(logtab))
    tau_res$Tau <- apply(logtab, 1, function(row) calculateTau(row))
    hist(tau_res$Tau, col="orange", main=paste("Tau index at", i,",", age),
         xlab= 'Tau index')
  }
  par(mfrow=c(1,1))
  dev.off()
}

## generating new histograms and venn diagrams on calculated Tau
## 
getNOhousekeeping <- function(age){
  tauokgenes <- c()
  myregex = paste0(age,".*\\.txt$")
  for (i in list.files("Tau/", pattern=myregex)){
    itext <- read.table(paste0("Tau/",i), sep='\t', header=T, row.names = 1)
    tmpvec <- itext %>% dplyr::filter(class !='housekeeping') %>% pull(id)
    tauokgenes <- c(tauokgenes, tmpvec)
  }
  return(unique(tauokgenes))
}


if (ulterior){
  tauokgenes.old <- getNOhousekeeping("Old")
  tauokgenes.young <- getNOhousekeeping("Young")
  venntau <- venn.diagram(x=list("OLD"=tauokgenes.old,
                                 "YOUNG"=tauokgenes.young),
                          height=15, width=10, units="cm",
                          col="transparent", fontfamily = "sans",
                          fill=c("darkgray", "#E69F00"),
                          alpha = 0.5, cex=.8, 
                          cat.col = c("black", "#E69F00"),
                          cat.cex = .8, margin=0.1,
                          print.mode=c("raw","percent"), filename=NULL)
  
  pdf(paste0(resdir,"tau_shared.pdf"),width=5, height=4)
  grid.arrange(gTree(children=venntau), 
               top=paste0("Specific or intermediate genes \n (by Tau >= 0.5)"),
               padding=unit(0.2, "line"))
  dev.off()
  #days = c("D0","D2","D4","D7")
  days = c("D4")
  day = "D4"
  for (day in days){
    oldtab <- read.table(paste0("Tau/TauSpecificity_",'Old',day,".txt"), sep='\t',
                       header = T)
    youngtab <- read.table(paste0("Tau/TauSpecificity_",'Young',day,".txt"), sep='\t',
                         header = T)
    # use ?scale_colour_colorblind
  }
}
###  For each cell type and days post injury, what is the distribution of specific and intermediate genes between young and old condition expression
getNOhousekeeping2 <- function(age,typeCell,Day){
  tauokgenes <- c()
  itext <- read.table(paste0("Tau/TauSpecificity_",age,Day,".txt"), sep='\t', header=T, row.names = 1)
  geneAgeDayCellSpe <- itext %>% dplyr::filter(class !='housekeeping') %>% dplyr::filter(str_detect(whichMAX,typeCell))
  tmpvec <- geneAgeDayCellSpe %>% pull(id)
  tauokgenes <- c(tauokgenes, tmpvec)
  return(unique(tauokgenes))
}

cellcolors = list(
  "ECs"="#0072B2",
  "FAPs"="#F0E442",
  "M1" = "#D55E00",
  "M2" =  "#CC79A7",
  "Neutro" =  "#009E73",
  "sCs" = "#56B4E9" )

tauokgenes.old <- getNOhousekeeping2("Old","ECs","D0")
tauokgenes.young <- getNOhousekeeping2("Young","ECs","D0")

DiagrammVennFacet <- function (typeCell,Day){
  tauokgenes.old <- getNOhousekeeping2("Old",typeCell,Day)
  tauokgenes.young <- getNOhousekeeping2("Young",typeCell,Day)
  #tauokgenes.old <- getNOhousekeeping2("Old","ECs","D0")
  #tauokgenes.young <- getNOhousekeeping2("Young","ECs","D0")
  col<- cellcolors[typeCell][[1]]
  venntau <- venn.diagram(x=list("OLD"=tauokgenes.old,
                                 "YOUNG"=tauokgenes.young),
                          height=10, width=5, units="cm",
                          col="transparent", fontfamily = "sans",
                          fill=c("darkgray", col),
                          alpha = 0.5, cex=.8, 
                          cat.col = c("black", col),
                          cat.cex = .8, margin=0.1,
                          print.mode=c("raw","percent"), filename=NULL
                          #, main =paste0(typeCell,Day," Specific or intermediate genes \n (by Tau >= 0.5)") 
                          )

  return(venntau)
}

pdf(paste0(resdir,"Old_Youg_DiagVenn_byTypeCellDays.pdf"))
par(mfrow = c(6,4) )
celltypev=c("ECs","FAPs","sCs")
Daysv=c("D0","D2","D4","D7")
DVallDays<-list()
for (typeCell in celltypev){
  for ( Day in Daysv) {
    CellDaysVennDiagramm<- DiagrammVennFacet(typeCell,Day)
    DVallDays<-list(DVallDays,CellDaysVennDiagramm)
    grid.arrange(gTree(children=CellDaysVennDiagramm))
  }
}

celltypev=c("M2")
Daysv=c("D2","D4","D7")
DVM2<-list()
for (typeCell in celltypev){
  for ( Day in Daysv) {
    CellDaysVennDiagramm<- DiagrammVennFacet(typeCell,Day)
    DVM2<-list(DVM2,CellDaysVennDiagramm)
    grid.arrange(gTree(children=CellDaysVennDiagramm))
  }
}
celltypev=c("M1")
Daysv=c("D2","D4")
DVM1<-list()
for (typeCell in celltypev){
  for ( Day in Daysv) {
    CellDaysVennDiagramm<- DiagrammVennFacet(typeCell,Day)
    grid.arrange(gTree(children=CellDaysVennDiagramm))
    DVM1<-list(DVM1,CellDaysVennDiagramm)
  }
}
CellDaysVennDiagramm<- DiagrammVennFacet("Neutro","D2")
DVNeutro<-list(CellDaysVennDiagramm)
grid.arrange(gTree(children=CellDaysVennDiagramm))
dev.off()
blank <- grid.rect(gp=gpar(col="white"))
facetPlot<-grid.arrange(gTree(children=DVallDays[[1]][[1]][[1]][[1]][[1]][[1]][[1]][[1]][[1]][[1]][[1]][[2]])
            ,gTree(children=DVallDays[[1]][[1]][[1]][[1]][[1]][[1]][[1]][[1]][[1]][[1]][[2]])
            ,gTree(children=DVallDays[[1]][[1]][[1]][[1]][[1]][[1]][[1]][[1]][[1]][[2]])
            ,gTree(children=DVallDays[[1]][[1]][[1]][[1]][[1]][[1]][[1]][[1]][[2]])
            ,gTree(children=DVallDays[[1]][[1]][[1]][[1]][[1]][[1]][[1]][[2]])
            ,gTree(children=DVallDays[[1]][[1]][[1]][[1]][[1]][[1]][[2]])
            ,gTree(children=DVallDays[[1]][[1]][[1]][[1]][[1]][[2]])
            ,gTree(children=DVallDays[[1]][[1]][[1]][[1]][[2]])
            ,gTree(children=DVallDays[[1]][[1]][[1]][[2]])
            ,gTree(children=DVallDays[[1]][[1]][[2]]),
            gTree(children=DVallDays[[1]][[2]]),
            gTree(children=DVallDays[[2]]),
            blank,
            gTree(children=DVM2[[1]][[1]][[2]]),
            gTree(children=DVM2[[1]][[2]]),
            gTree(children=DVM2[[2]]),
            blank,
            gTree(children=DVM1[[1]][[2]]),
            gTree(children=DVM1[[2]]),
            blank,
            blank,
            gTree(children=DVNeutro[[1]]),
            blank,
            blank,
            ncol=4,nrow=6)
save_plot("Tau/Old_Youg_DiagVenn_byTypeCellDay_nolegend.png",facetPlot,base_height = 15,base_width = 12)


## autre représentation de la spécificité des gènes entre Y et O
vectTypeCell<-c("ECs","FAPs","M1","M2","Neutro","sCs")
vectDay<-c("D0","D2","D4","D7")
deesignTypeCellDays<-data.frame("D0"=c(TRUE,TRUE,FALSE,FALSE,FALSE,TRUE),"D2"=c(TRUE,TRUE,TRUE,TRUE,TRUE,TRUE),"D4"=c(TRUE,TRUE,TRUE,TRUE,FALSE,TRUE),"D7"=c(TRUE,TRUE,FALSE,TRUE,FALSE,TRUE))
row.names(deesignTypeCellDays)<-vectTypeCell
countIntersectYO <- tibble(NumberGeneSpe=numeric(),PercentGeneSpe=numeric(),DaypostInjury=character(),TypeCell=character(),Group=character())
for (typeCell in vectTypeCell){
  print(typeCell)
  for (Day in vectDay){
    if ( isTRUE(deesignTypeCellDays[typeCell,Day])){
      print(Day)
      tauokgenes.old <- getNOhousekeeping2("Old",typeCell,Day)
      tauokgenes.young <- getNOhousekeeping2("Young",typeCell,Day)
      totuniqgene<-length(unique(c(tauokgenes.old,tauokgenes.young)))
      countIntersectYO<- countIntersectYO %>% add_row(NumberGeneSpe=length(setdiff(tauokgenes.young,tauokgenes.old)),PercentGeneSpe=length(setdiff(tauokgenes.young,tauokgenes.old))/totuniqgene*100,DaypostInjury=Day,TypeCell=typeCell,Group="OnlyYoung")
      countIntersectYO<- countIntersectYO %>% add_row(NumberGeneSpe=length(setdiff(tauokgenes.old,tauokgenes.young)),PercentGeneSpe=length(setdiff(tauokgenes.old,tauokgenes.young))/totuniqgene*100,DaypostInjury=Day,TypeCell=typeCell,Group="OnlyOld")
      countIntersectYO<- countIntersectYO %>% add_row(NumberGeneSpe=length(intersect(tauokgenes.young,tauokgenes.old)),PercentGeneSpe=length(intersect(tauokgenes.young,tauokgenes.old))/totuniqgene*100,DaypostInjury=Day,TypeCell=typeCell,Group="ShareYoungOld")
      
    }
  }
}
countIntersectYO<-countIntersectYO %>% arrange(Group,DaypostInjury)
typebytime<-ggplot(countIntersectYO, aes(x=DaypostInjury,y=PercentGeneSpe,fill=Group,color=Group,group=Group)) +
  geom_point(alpha=0.8)+
  geom_line()+
  scale_color_manual(values=c("#999999", "#E69F00", "#56B4E9"))+
  theme_ipsum()+
  facet_grid(vars(TypeCell)) +
  ggtitle("Gene Specificity Share for Young and Old by typeCell")

typebytime
save_plot("Tau/Old_Youg_percentgenespecificity_byTypeCellDay
          .png",typebytime,base_height =10,base_width = 8)
