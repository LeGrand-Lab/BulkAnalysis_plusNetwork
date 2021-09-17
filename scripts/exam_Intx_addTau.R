# test to cross specificity and DE results
# by day by cell type
# JohaGL
setwd("~/BulkAnalysis_plusNetwork/")
library(tidyverse)
odir = "exam_INTER_conditions/static/"
DEf = readRDS(paste0(odir,"shot_rds_full.rds"))
dim(DEf)
head(DEf,2)

toyds <- sample_n(DEf, 2000)
consensus_tau <- list()
daysl = c("D0","D2","D4","D7")
for (day in daysl){
  ty <- read.table(paste0("Tau/TauSpecificity_Young",day, ".txt"), sep="\t", header=T)
  to <- read.table(paste0("Tau/TauSpecificity_Old",day, ".txt"), sep="\t", header=T)
  # remember = Tau is calculated on meanTPM (mean over replicates' TPMs)
  #  therefore, its value is influenced by changes in expression due to age
  # for example, a hypothetical gene1 : [FAPS, sCs, ECs, M1], their meanTPMs
  # being :   Old =  [100,100,100,100]  ; Young = [100,2000,100,100]
  #  the result will say that gene1 is 'sCs' specific, but only in Young.
  #  As a matter of simplicity, by day and by cell type:
  #             if gene has Tau > 0.5 at least in one of the ages
  # it will be included in consensus list:
  brut <- rbind(ty,to)
  brut <- brut  %>% filter(Tau > 0.5) %>% filter(nbMAX == 1 ) %>%
  group_by(symbol, whichMAX) %>%    # whichMAX is the cell type
  slice(which.max(Tau)) 
  consensus_tau[[day]] <- brut %>% select(symbol, Tau, class, whichMAX)
}
for (i in (names(consensus_tau))){
  print("")
  print(i)
  print(table(consensus_tau[[i]]$whichMAX))
  print(table(consensus_tau[[i]]$class))
  print("-----------------")
}

saveRDS(consensus_tau, paste0(odir,"conseTau4DEGs/", "conseTau.rds"))
# 
# [1] ""
# [1] "D0"
# 
# ECs FAPs  sCs 
# 1265 2755 1682 
# 
# intermediate     specific 
# 3020         2682 
# [1] "-----------------"
# [1] ""
# [1] "D2"
# 
# ECs   FAPs     M1     M2 Neutro    sCs 
# 1463   2417    454   1359   1906   2053 
# 
# intermediate     specific 
# 4914         4738 
# [1] "-----------------"
# [1] ""
# [1] "D4"
# 
# ECs FAPs   M1   M2  sCs 
# 1573 3104  795 1008 1559 
# 
# intermediate     specific 
# 3990         4049 
# [1] "-----------------"
# [1] ""
# [1] "D7"
# 
# ECs FAPs   M2  sCs 
# 1219 3379 1367 1887 
# 
# intermediate     specific 
# 3841         4011 
# [1] "-----------------"

