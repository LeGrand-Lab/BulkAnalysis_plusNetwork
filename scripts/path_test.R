BiocManager::install("pathview")
library("pathview")
data(gse16873.d)
data(demo.paths)
data(paths.hsa)

i <- 1
pv.out <- pathview(gene.data = gse16873.d[, 1], pathway.id =
                     demo.paths$sel.paths[i], species = "hsa", out.suffix = "gse16873",
                   kegg.native = TRUE)

myvec <- localgenes$logFC
names(myvec) <- localgenes$id
# note 'Metabolic pathways' is huge!!! no visually informative
code = rownames(tmpkegg_t_t[tmpkegg_t_t$Pathway=="NOD-like receptor signaling pathway",] )

pv.out <- pathview(
  gene.data = myvec,
  pathway.id = str_replace(code,"path:",""),
  species = "mmu",
  out.suffix = str_replace(code,"path:",""),
  gene.idtype="entrez"
)
