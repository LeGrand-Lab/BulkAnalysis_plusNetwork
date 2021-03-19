####
setwd("~/bulk_analysis")
listgenes <- readRDS("transf_files/listGenes_DEG_by_cluster_by_celltype.Rds")

# replace id_ensemble by symbols in matrices rownames...
ids = rownames(fmat)
gdf = read.table("genesinfo.csv")
symbols = c()
for (i in 1:length(ids)){
  symbols <- c(symbols, gdf[gdf$Geneid==ids[i],]$symbol)
}
symbols <- make.unique(symbols)
(c(head(symbols),"...",tail(symbols)))
rownames(fmat) <- symbols
rownames(fTPM) <- symbols

# ####
genes.under <-  names(meanTPM[meanTPM<=0.01])
genes.over.cutoff = names(meanTPM[meanTPM>0.01])

sup_m <- filtered.mat[genes.over.cutoff,]
sup_t <- filteredTPM[genes.over.cutoff,]

k_m <- apply(sup_m, 1, sum)
k_m[k_m<=1]

k2_m <- apply(sup_m, 1, sum)

k_m[k_m<=1]

inf_t <- filteredTPM[genes.under,]
inf_m <- filtered.mat[genes.under,]

i2_m <- apply(inf_m, 1, function(x) if(sum(x>10)==127){return(TRUE)}
              else{return(FALSE)})
i2_m[i2_m==T]

# ENSMUSG00000064341 ENSMUSG00000064345 ENSMUSG00000064363 
# TRUE               TRUE               TRUE 
