##
# iTALK on D7 data, old vs young
##
# iTALK requires
# $ sudo apt-get install libcairo2-dev
# $ sudo apt-get install libnode-dev
# devtools::install_github("Coolgenome/iTALK", build_vignettes = TRUE)
library(dplyr)
library(tidyverse)
library(dplyr)
library(tidyverse)
library(ggplot2)
library(cowplot)
library(forcats)
library(RColorBrewer)
library(DESeq2)
library(gridExtra)
library(reshape2)
library(ggthemes)
library(iTALK)
# demo by iTALK
setwd("~/bulk_analysis/")
odir = "italk/"
symbols <- 'data/genesinfo.csv'
prefil_cou <- "data/prefiltered_counts.rds"
prefil_tpm <- "data/prefiltered_TPM.rds"
metadata.rds <- "data/metadata.rds"

## import dataframes / matrices
symb_df <- read.table(symbols, header=T)
fmat <- readRDS(prefil_cou)
metadata <- readRDS(metadata.rds)
metadata <- metadata %>% mutate(timetype=paste0(time,".",type)) 
## define day !!: 
day = "D7"
## get subset:
xmat <- fmat[,str_detect(colnames(fmat),day)]
xmeta <- metadata[str_detect(rownames(metadata),day),]
x.keep <- apply(xmat, 1, function(row) ifelse(sum(row >=5)>= 3, TRUE, FALSE))
xmat <- xmat[x.keep,]

xmat
