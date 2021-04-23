# plot genes into pathview style
# johaGL 2021
#BiocManager::install("pathview")
library("pathview")
data(gse16873.d)
data(demo.paths)
data(paths.hsa)
setwd("~/BulkAnalysis_plusNetwork")
resdir = "path_classic_tools/"

vrs = "vB"
age="Young"
### load Dyn edger genes AND kegg info
DYNgenes <- read.table(paste0(resdir,"edger_dynINTRA_",age,vrs,".txt"), sep='\t',
                       header=T)
gokegg <- readRDS( file=paste0(resdir,"edger_dynINTRAkeggo",age,vrs,".rds"))
# do match both infos
# pick celltype and contrast
ct = "ECs"
localgenes = DYNgenes %>% filter(type==ct) 
names(gokegg[[paste0(ct,"_kegg")]]) # the contrasts tested : "D2vsD0" "D4vsD2" "D7vsD4"
contrast = "D7vsD4"
head(gokegg[[paste0(ct,"_kegg")]][[contrast]])
tmpkegg_t_t = gokegg[[paste0(ct,"_kegg")]][[contrast]]
dim(tmpkegg_t_t)

###
myvec <- localgenes$logFC
names(myvec) <- localgenes$id
code = rownames(tmpkegg_t_t[tmpkegg_t_t$Pathway=="NOD-like receptor signaling pathway",] )

dir.create(paste0(resdir,age,ct,contrast),recursive=T)
setwd(paste0(resdir,age,ct,contrast))
pathview(
  gene.data = myvec,
  pathway.id = str_replace(code,"path:",""),
  species = "mmu",
  out.suffix = str_replace(code,"path:",""),
  gene.idtype="entrez"
)
setwd("~/BulkAnalysis_plusNetwork")
