### demo italk
# iTALK requires
# $ sudo apt-get install libcairo2-dev
# $ sudo apt-get install libnode-dev
library(dplyr)
library(tidyverse)
# devtools::install_github("Coolgenome/iTALK", build_vignettes = TRUE)
library(iTALK)
# demo by iTALK
setwd("~/bulk_analysis/")
demodata = read.table("italk_example_data.txt", sep='\t', header=T,
                      stringsAsFactors = F)
highly_exprs_genes<-rawParse(demodata,top_genes=50,stats='mean')
# find the ligand-receptor pairs from highly expressed genes
comm_list<-c('growth factor','other','cytokine','checkpoint')
cell_col<-structure(c('#4a84ad','#4a1dc6','#e874bf','#b79eed', '#ff636b', '#52c63b','#9ef49a'),
                    names=unique(demodata$cell_type))
par(mfrow=c(1,2))
res<-NULL
for(comm_type in comm_list){
  res_cat<-FindLR(highly_exprs_genes,datatype='mean count',comm_type=comm_type)
  res_cat<-res_cat[order(res_cat$cell_from_mean_exprs*res_cat$cell_to_mean_exprs,decreasing=T),]
  #plot by ligand category
  #overall network plot
  NetView(res_cat,col=cell_col,vertex.label.cex=1,arrow.width=1,edge.max.width=5)
  #top 20 ligand-receptor pairs
  LRPlot(res_cat[1:20,],datatype='mean count',cell_col=cell_col,
         link.arr.lwd=res_cat$cell_from_mean_exprs[1:20],
         link.arr.width=res_cat$cell_to_mean_exprs[1:20])
  title(comm_type)
  res<-rbind(res,res_cat)
}
res<-res[order(res$cell_from_mean_exprs*res$cell_to_mean_exprs,decreasing=T),][1:20,]
NetView(res,col=cell_col,vertex.label.cex=1,arrow.width=1,edge.max.width=5)
LRPlot(res[1:20,],datatype='mean count',cell_col=cell_col,
       link.arr.lwd=res$cell_from_mean_exprs[1:20],
       link.arr.width=res$cell_to_mean_exprs[1:20])

# do diff exp analysis
data <- demodata
# randomly assign the compare group to each sample
data <- data %>% mutate(compare_group=sample(2,nrow(data),replace=TRUE))
deg_nk <- DEG(data,
            method='MAST',contrast=c(2,1))
data %>% group_by(cell_type,compare_group) %>% summarise(pop=n())
tmp = data %>% dplyr::select("cell_type","compare_group")


## run DEG using deseq2 in simulated counts:
simcounts <- as.data.frame(matrix(sample.int(700,size = 8*50, replace=T),
                                  nrow=8,ncol=50))
colnames(simcounts) <- as.character(seq(1,50))
simcounts$cell_type <- c("myo","myo","myo", "myo", "wtf","wtf","wtf" ,"wtf")
# here contrast 'wtf' against all others
simcounts$compare_group  <- factor(c(1,1,2,2,1,1,2,2))

# test 1 vs 2 (2 is reference )
for (cellty in unique(simcounts$cell_type)){
  deg_test <- DEG( simcounts %>% filter(cell_type==cellty),
                   method='DESeq2',contrast=c(1,2))
  print(head(deg_test))
}

# do for each celltype:

