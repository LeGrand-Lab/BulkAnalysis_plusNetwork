## first make sure that required system libs, run in Debian, Ubuntu, etc:
# apt-get install libssl-dev 
# apt-get install libcurl4-openssl-dev
# apt-get install libxml2-dev

print("you can say 'no' to question : Do you want to install from sources the packages which need compilation?" )

# ##
# Installing CRAN packages
# ##
cranpkg <- c("Rcpp", "curl", "xml2", "stringr", "dplyr","tidyverse", "ggplot2", 
	"Matrix", "reshape2", "gridExtra", "matrixStats",
	"boot", "class", "nnet", #required for DropletUtils
	"future", "future.apply", "glmpca", 
	"gprofiler2", "msigdbr",
	"MASS", "sctransform", "Seurat",
	"reticulate", # for python interface
	"rmarkdown", 
	"igraph", "shiny", "networkD3", "kml")

# installing cran  packages
sapply(cranpkg, function(x) {
		if (!requireNamespace(x, quietly=F)) {
			install.packages(x)
		} 
  })

# check everything installed correctly: returns TRUE if each package in the list is already installed, FALSE otherwise
resultsinstall1 <- lapply(cranpkg, require, character.only=TRUE, quietly=F)

print("printing CRAN packages' names that could not be installed:")
notinstalledcran <- c()
names(resultsinstall1) <- cranpkg
for (k in names(resultsinstall1)){
  if (!resultsinstall1[[k]] ){
    print(paste(k, " : could not be installed"))
    notinstalledcran  <- c(notinstalledcran, k)
  }
}


# ##
# installing Biocmanager or devtool packages:
# ##


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
devtools::install_github("vqv/ggbipplot")

bioc_or_dev <- c("AnnotationDbi","org.Mm.eg.db", "harmony", "SingleCellExperiment", 
                 "scater", , "DESeq2", 
                 "scran", "DropletUtils", "clusterProfiler",
                 "gage", "fgsea", "forcats", "UpSetR","ggpubr","kml", "ggplotify",
                 "ComplexHeatmap")
	
installbd(bioc_or_dev) # use function to install
# verify that all bio or devtools installed correctly: 
resultsinstall2 <- lapply(bioc_or_dev, require, character.only=TRUE, quietly=F)
print("printing Bioconductor or devtools packages that could not be installed")
notinstalledbd <- c()
names(resultsinstall2) <- bioc_or_dev
for (k in names(resultsinstall2)){
  if (!resultsinstall2[[k]] ){
    print(paste(k, ": could not be installed"))
    notinstalledbd <- c(notinstalledbd, k)
  }
}


# sapply(bioc_or_dev, function(x){
#   tryCatch({install.packages(x)
#     return(paste(x," ok "))},
#     error=function(cond){
#       message(cond)
#       message(" ** trying with bioconductor ** ")
#       BiocManager::install(x)
#     })
# })

## END
# note: a package installed from unknown source : "HSMMSingleCell"
