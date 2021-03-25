## first make sure that required system libs, run in Debian, Ubuntu, etc:
# apt-get install libssl-dev 
# apt-get install libcurl4-openssl-dev
# apt-get install libxml2-dev

cranpkg <- c("Rcpp", "curl", "xml2", "stringr", "dplyr","tidyverse", "ggplot2", 
	"Matrix", "reshape2", "gridExtra", "matrixStats",
	"boot", "class", "nnet", #required for DropletUtils
	"future", "future.apply", 
	"MASS", "sctransform", "Seurat",
	"reticulate", # for python interface
	"rmarkdown", 
	"igraph", "shiny", "networkD3")

# installing cran or other R packages
sapply(cranpkg, function(x) {
		if (!requireNamespace(x, quietly=F)) {
			install.packages(x)
		} 
  })

# check everything installed correctly
lapply(cranpkg, require, character.only=TRUE, quietly=F)

# installing Biocmanager packages
install.packages("BiocManager")
install.packages("devtools")
remotes::install_github('chris-mcginnis-ucsf/DoubletFinder')
installbd <- function(pks){
	sapply(pks,
		if (!requireNamespace(pks, quietly=F))
		tryCatch({
			BiocManager::install(pks)
		},error = function(e) {
			devtools::install_github(pks)
		})
	)	
  }	

bioc_or_dev <- c("harmony", "SingleCellExperiment", "scater", 
                 "scran", "DropletUtils", "clusterProfiler")
	
installbd(bioc_or_dev)
lapply(bioc_or_dev, require, character.only=TRUE, quietly=F)

libs4bulk <- c("forcats", "UpSetR","ggpubr","kml", "ggplotify")

sapply(libs4bulk, function(x){
  tryCatch({install.packages(x)
    return(paste(x," ok "))},
    error=function(cond){
      message(cond)
      message(" ** trying with bioconductor ** ")
      BiocManager::install(x)
    })
})

# note: a package installed from unknown source : "HSMMSingleCell"
