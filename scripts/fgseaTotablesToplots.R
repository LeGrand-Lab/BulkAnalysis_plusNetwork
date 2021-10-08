# Organize GSEA results into tables
#   check fgsea_dayandcelltype.R
# output dir: exam_INTER_conditions/static/GSEA/ 
# and pdf plots
# --
# johaGL
library(tidyverse)
library(ComplexHeatmap)
library(msigdbr)
library(cowplot)

setwd("~/BulkAnalysis_plusNetwork/")
odir <- "exam_INTER_conditions/static/"
tablesneeded <- F
mat4heatmapneeded <- F
gseaoutfull = "fgsea_bd_bct_full.rds" 
gseaoutfiltered = "fgsea_bd_bct_filtered.rds"

fullrds <- readRDS(paste0(odir,"GSEA/rds/", gseaoutfull))
pathsFiltered = readRDS( paste0(odir, "GSEA/rds/",  gseaoutfiltered) )

thegmt <- msigdbr(species = "Mus musculus", 
                  category = 'C2', 
                  subcategory=c('CP:REACTOME'))
msigdbr_list = split(x = thegmt$gene_symbol, f = thegmt$gs_name)

fullDEsta = readRDS(paste0(odir, "rds/shot_rds_full.rds"))

if (tablesneeded){
  system(paste0("cd ",odir,"GSEA/; ", 
               "if [ ! -d csv/ ]; then mkdir csv; else echo 'csv/ exists, no action';fi"))
  print(" ** see keys in full results rds list of lists ** ")
  print(names(fullrds)) # the days
  print(names(fullrds[["D0"]])) # the cell types (that exist in given day)
  days = names(fullrds)
  for (d in names(fullrds)){
    print(paste("\n^^^^^", d, "^^^^^"))
    for (t in names(fullrds[[d]])){
      print(paste("saving", d, t, "as csv"))
      print(t)
      print(dim(fullrds[[d]][[t]]))
      tmpdf <- fullrds[[d]][[t]]
      
      tmpdf <- tmpdf %>% 
        mutate(gene_symbols = sapply(leadingEdge, function(x) paste(x, collapse=", ")),
               .before = leadingEdge) %>%
        select(-leadingEdge)
      
      write.table(tmpdf, paste0(odir,"GSEA/csv/", t,"_",d,"_pathways.csv"), sep="\t", 
                  col.names = T, row.names=F              )
    }
  }
  rm(fullrds)
}

## NOTE : tryCatch in innerplot_fun allowed to handle the error : 
#  Error in calcGseaStat(statsAdj, selectedStats = pathway, returnAllExtremes = TRUE) : 
# GSEA statistic is not defined when all genes are selected 
 
innerplot_fun <- function(fgsea.dayhere, d, CT, NumP, msigdbr_list){
  print("running innerplot_fun")
  fgseaRes <- fgsea.dayhere[[CT]] %>% mutate(log2err=replace_na(log2err, 1))  ## not use drop_na : all down faps D2 are log2err NA
  print(min(fgseaRes$padj))
  topPathwaysUp <- fgseaRes[ES>0][head(order(padj),n=NumP), ]
  print(paste("nb of upregul paths", length(topPathwaysUp)))
  topPathwaysDown <- fgseaRes[ES<0][head(order(padj), n=NumP), ]
  print(paste("nb of downreg paths", length(topPathwaysDown)))
  topPathBoth <- c(topPathwaysUp$pathway, rev(topPathwaysDown$pathway))
  a <- sort(getgeneslist_mod(topPathwaysUp, d), decreasing = T)
  b <- sort(getgeneslist_mod(topPathwaysDown, d), decreasing = T)
  gseagenes <- c(a, b)
  gseagenes <- sort(gseagenes, decreasing=T)
  print(gseagenes)
  ouif = fgsea::plotGseaTable(msigdbr_list[topPathBoth], gseagenes, fgseaRes, 
            gseaParam = 0.5 ,  colwidths = c(3, 3, 0.8, 1.2, 1.2),  render=F) 
  print('doneplotgseatable')
  plotsenrichu_ <- list()
  for (i in 1:NumP){
    pup = topPathwaysUp[i,]
    print("a path up")
    print(pup)
    tmpup <- NULL
    tryCatch({ 
      tmpup <- fgsea::plotEnrichment(msigdbr_list[[pup$pathway]], gseagenes) + 
      labs(title = str_replace(pup$pathway, "REACTOME_", ""),
           caption=paste( "(NES:", round(pup$NES, 2), ", padj :", 
                          round(pup$padj,2),")") )
    }, 
    warning = function(w) {}, error = function(e){print("paths up error")} ) # endtrycatch
    if ((dim(tmpup$data)[1] <= 2) || is.null(tmpup) ){
      plotsenrichu_[[i]] <- NULL
    }else{ plotsenrichu_[[i]] <- tmpup + theme(plot.title = element_text(size=8))}
  }
  print("done enrichu")
  plotsenrichdw_ <- list()
  for (i in 1:NumP){
    pdw = topPathwaysDown[i,]
    print("b path down")
    print(pdw)
    tmpdw <- NULL
    tryCatch({
      tmpdw <- fgsea::plotEnrichment(msigdbr_list[[pdw$pathway]], stats=gseagenes) + 
          labs(title = str_replace(pdw$pathway,"REACTOME_", ""),
           caption=paste( "(NES:", round(pdw$NES,2), ", padj :", 
                          round(pdw$padj,2),")") )
        }, 
      warning = function(w){}, error = function(e){print("paths down error")} ) # trycatch
    if (dim(tmpdw$data)[1] <= 2 || is.null(tmpdw)){
      plotsenrichdw_[[i]] <- NULL
    } else {plotsenrichdw_[[i]] <- tmpdw + theme(plot.title = element_text(size=8)) }
  }
  print("done down")
  print(
  plot_grid(
    # inner plotgrid 1 : the 'plotGseaTable',
    plot_grid(ggdraw() + draw_label(paste(d, CT, ": Top enriched Pathways (GSEA), Old vs Young")),
              plot_grid(NULL, ouif, NULL, nrow=1, rel_widths = c(10,9,1)), # NULL elem helps de-truncate names
              nrow=2, rel_heights = c(2, 10)) ,
    # inner plotgrid 2 : the 'plotEnrichment' : u_ == UP, dw_ == DOWN
    plot_grid(
      plot_grid(plotlist = plotsenrichu_, nrow = 1, ncol=NumP+1),
      plot_grid(plotlist = plotsenrichdw_,  nrow = 1, ncol=NumP+1) , nrow = 2, rel_heights = c(5,5)), 
    
    nrow=2, rel_heights = c(4,5)
    )# end plotgrid all  
  ) # end print
  print("done plotgrid")
}

getgeneslist_mod <- function(gseadatafr, d){
  outi <- c()
  for (i in gseadatafr$leadingEdge){
    outi <- c(outi, i)
  }
  moo <- unique(outi)
  tmpdfdeg <- fullDEsta %>% filter(day==d)
  lfcs <- tmpdfdeg[match(moo, tmpdfdeg$symbol),]$log2FoldChange
  names(lfcs) <- moo
  return(lfcs[!is.na(lfcs)])
}

plotme_mod <- function(pathsFiltered, d, outfileprefix, NumP=5){
  fgsea.dayhere <- pathsFiltered[[d]]
  print(paste("pathways object is a list of lists, inner list has celltypes:",
        ifelse((typeof(fgsea.dayhere) == "list" & 
                 is.null(colnames(fgsea.dayhere))), "ok",  #colnames(a list) is null
                  "ERROR no celltypes in names(fgsea.dayhere), type prblm" )))
  types = names(fgsea.dayhere)
  print(types)
  pdf(paste0(odir, "GSEA/", outfileprefix, d,".pdf"), width = 14, height=8)
  for (CT in types){
    print(CT)
    print( innerplot_fun( fgsea.dayhere, d, CT, NumP, msigdbr_list) )
  }  
  dev.off()
}
rm(d)
plotme_mod(pathsFiltered, 'D0', "bd_bct_result_", NumP=5); dev.off()
plotme_mod(pathsFiltered, 'D2', "bd_bct_result_", NumP=5); dev.off()
plotme_mod(pathsFiltered, 'D4', "bd_bct_result_", NumP=5) ; dev.off()
plotme_mod(pathsFiltered, 'D7', "bd_bct_result_", NumP=5) ; dev.off()



# ============= Prepare matrices for individual (dummy) plots ==================
## prepare matrices for heatmaps (small heatmaps )
if (mat4heatmapneeded){
  m_ = list()
  for (k in c('D0','D2', 'D4', 'D7')){
    cts <- unique(names(pathsFiltered[[k]])) 
    m_[[k]] <- list()
    for (CT in cts){
      m_[[k]][[CT]] <- list("UP"=NULL,"DOWN"=NULL)
      for(S in c("UP","DOWN")){
        here.data = as.data.frame(pathsFiltered[[k]][[CT]])
        here.data <- here.data %>% filter(sens == S) %>% 
          mutate(path4graph = str_replace(pathway, "REACTOME_", "")) 
        here.dico = list()
        if(length(here.data$path4graph) == length(unique(here.data$path4graph) )){
          print(paste("ok, processing", k, CT, S))
          rownames(here.data) =  here.data$path4graph
          aggreg_genes = c()
          for (path in here.data$path4graph){
            tmpstr <- str_replace_all(here.data[path,]$leadingEdge, 
                                      c('c\\("' = '',   '"' = '', '\\)' = ''))
            here.dico[[path]] <- unname(unlist(str_split(tmpstr,", ")))
            aggreg_genes <- c(aggreg_genes, here.dico[[path]] )
          }
          colnamesgenes = unique(aggreg_genes)
          matrx = array(NA, dim=c(length(here.data$path4graph), length(colnamesgenes)))
          rownames(matrx) = here.data$path4graph
          colnames(matrx) = colnamesgenes
          for (path in here.data$path4graph){
            localgenes <- here.dico[[path]]
            for (g in localgenes){
              matrx[path,g] <-  here.data[path,]$NES # if more than once, yields 2 or more
            }
          }
          # clear rows and columns containing 2 or less values
          matrx <- matrx[rownames(matrx) != "DISEASE",] # take away "DISEASE" ! 
          rtokeep <- apply(matrx, 1, function(x) sum(!is.na(x)) >= 3)
          newm <- matrx[rtokeep,]
          m_[[k]][[CT]][[S]] <- newm
        } else {
          print(paste("error, this", k, CT ,"  contains repeated pathwaynames") )
          stop()
        } # endif pathways vectors are unique (no elements repeated)
      } # end for sens UP or DOWN
    }# end for CT in cts
  }# end for k in vector of days
  saveRDS(m_, file=paste0(odir,"GSEA/rds/fgsea_matrices4heatmaps.rds" ))
} # end if mat4heatmapneeded

# ===================== Prepare giant single matrix ============================
#           and exclude DISEASE (use str_detect) and DEVELOPMENT 

givemeltedDataframe <- function(pathsFiltered, minnbg, nbminpadj, COLLECTION){
  r_path = c(); r_sens = c(); r_NES = c(); 
  r_gene = c(); r_celltype = c(); r_day = c(); r_padj = c()
  counter = 0
  for (d in names(pathsFiltered)){
    for (t in names(pathsFiltered[[d]])){
      tmp <- pathsFiltered[[d]][[t]] %>% filter(str_detect(pathway, COLLECTION)) %>%
          group_by(sens) %>% slice_min(padj, n=nbminpadj)
      path_ = tmp$pathway
      for (i in 1:length(path_)){
        pw = tmp[i,]$pathway
        NES = tmp[i,]$NES
        se = tmp[i,]$sens
        pj = tmp[i,]$padj
        if (str_detect(pw, "DISEASE") | str_detect(pw, "DEVELOPMENT") | str_detect(pw, "INFECTION") ){
          print("excluded disease, devel and infection terms")
        } else {
          for(lg in tmp[i,]$leadingEdge){
            if(length(lg) >= minnbg ){
              for (g in lg){
                print(g)
                counter = counter + 1
                r_path = c(r_path, pw)
                r_sens = c(r_sens, se) 
                r_NES = c(r_NES, NES)
                r_gene =  c(r_gene, g)
                r_celltype = c(r_celltype, t)
                r_day = c(r_day, d)
                r_padj = c(r_padj, pj)
              } 
            }
          }
        }
      }#end for i in 1:length(path_)
    }
  }
  dfresu <- data.frame("pathway" = r_path, "sens"= r_sens, "NES" = r_NES,
                            "gene" = r_gene, "celltype" = r_celltype,
                            "day"=r_day, "padj" = r_padj)
  dfresu <- dfresu %>% mutate(unigene = paste0(gene,"_", celltype),
                                        unipath = paste0(pathway,"_", day)) 

  return(dfresu)
} 

MELTED <- givemeltedDataframe(pathsFiltered, 5, 15, "REACTOME") 
names(MELTED) 
head(MELTED)
# path sens        NES   gene celltype day    unigene                                            unipath
# 1 REACTOME_ABC_FAMILY_PROTEINS_MEDIATED_TRANSPORT   UP  0.7683668  Abcb4      ECs  D0  Abcb4_ECs REACTOME_ABC_FAMILY_PROTEINS_MEDIATED_TRANSPORT_D0
# 2 REACTOME_ABC_FAMILY_PROTEINS_MEDIATED_TRANSPORT   UP  0.7683668  Abcc1      ECs  D0  Abcc1_ECs REACTOME_ABC_FAMILY_PROTEINS_MEDIATED_TRANSPORT_D0
# 3 REACTOME_ABC_FAMILY_PROTEINS_MEDIATED_TRANSPORT   UP  0.7683668  Abcc9      ECs  D0  Abcc9_ECs REACTOME_ABC_FAMILY_PROTEINS_MEDIATED_TRANSPORT_D0
# 4 REACTOME_ABC_FAMILY_PROTEINS_MEDIATED_TRANSPORT   UP  0.7683668  Abca5      ECs  D0  Abca5_ECs REACTOME_ABC_FAMILY_PROTEINS_MEDIATED_TRANSPORT_D0
# 5                 REACTOME_ADAPTIVE_IMMUNE_SYSTEM DOWN -1.0629117 Zbtb16      ECs  D0 Zbtb16_ECs                 REACTOME_ADAPTIVE_IMMUNE_SYSTEM_D0

print("adding information about genes (padj, lfc), from 'rds/shot_rds_full.rds'")
toto <- MELTED %>% select(unigene, gene, celltype, day)
fuu <- fullDEsta %>% select(symbol, padj, log2FoldChange, type, day) %>% 
  mutate(gene = symbol, celltype = type, genepadj = padj)
tototo <- left_join(toto, fuu, by = c("day", "gene", "celltype" ) )
tototo <- tototo  %>% select(unigene, genepadj, log2FoldChange, day )
ttp = left_join(MELTED, tototo, by=c("unigene", "day"))  %>% unique() 
MELTED = ttp

rm(toto, fuu, tototo, ttp)

approve <- MELTED %>% group_by(celltype, day) %>% slice_min(genepadj, n=10) %>%
  filter(abs(log2FoldChange)>=0.3)
length(unique(approve$pathway))
length(unique(approve$unigene))

MELTED <- approve
# funciton to build megamatrix:
dfTomatrix <- function(dfx){
  mat <- array(NA,dim=c( length(unique(dfx$unipath)),
                          length(unique(dfx$unigene))
  ))
  rownames(mat) <- unique(dfx$unipath)
  colnames(mat) <- unique(dfx$unigene)
  for (p in dfx$unipath){
    for (g in dfx$unigene){
      tryCatch({
        val <- dfx[dfx$unipath == p & dfx$unigene == g,]$NES
      },
      warning = function(w){
      },
      error = function(e){
      })
      if (length(val)!= 0 ){
        mat[p,g] <- val
        print(val)
      }
    }
  }
  return(mat)
}

print("building megamatrix, take a cup of coffee (takes 15 minutes)")
m_reac <- dfTomatrix(MELTED)
print("clearing megamatrix")
keepro <- apply(m_reac, 1, function(r) sum(!is.na(r))>1)
m_reac2 <- m_reac[keepro,]

# =============================== REACTOME plot ================================ 
# render plot reactome

reac_daysvec <- sapply(rownames(m_reac2), function(x) {
  ltmp <- unlist(str_split(x, "_")) 
  return(ltmp[length(ltmp)])})
reac_celltyvec <- sapply(colnames(m_reac2),function(x){
  ltmp <- unlist(str_split(x, "_"))
  return(ltmp[length(ltmp)])
})
reac_newrownames <-  make.unique( sapply( names(reac_daysvec), function(x){  
  # "REACTOME_COLLAGEN_FORMATION_D7"
  ltmp <- unlist(str_split(str_replace(x,"REACTOME_",""), "_"))  
  mot <- paste0(ltmp[1:(length(ltmp)-1)], collapse="_")
  return(mot)  # "COLLAGEN_FORMATION
}) )
reac_newcolnames <- make.unique( sapply(names(reac_celltyvec), function(x){
  ltmp <- unlist(str_split(x, "_"))  # "Lamc1_sCs"
  return(ltmp[1]) # Lamc1
}) )
# assigne new names to all objects
rownames(m_reac2) <- reac_newrownames
colnames(m_reac2) <- reac_newcolnames
names(reac_daysvec) <- reac_newrownames
names(reac_celltyvec) <- reac_newcolnames

reac_splitdays <- data.frame(x=reac_daysvec)
rownames(reac_splitdays) <- rownames(m_reac2)
reac_splitcols <- data.frame(y=reac_celltyvec)
rownames(reac_splitcols) <- colnames(m_reac2)

ha = HeatmapAnnotation(foo = anno_block(gp = gpar(fill = c("darkblue", "orange",
                                                           "firebrick","violet",
                                                           "darkgreen","royalblue") )
)) #  labels = c("ECs","FAPs","M1","M2","Neutro","sCs") 
oh <- ComplexHeatmap::Heatmap(m_reac2, 
                              na_col = "whitesmoke",
                              cluster_rows = FALSE,
                              cluster_columns = FALSE,
                              row_split = reac_splitdays,
                              column_split = reac_splitcols,
                              row_names_gp = gpar(fontsize = 9),
                              column_names_gp = gpar(fontsize = 8),
                              name = "REACTOME",
                              column_names_rot = 45,
                              heatmap_legend_param = list(title="NES",
                                                          direction = "horizontal"),
                              top_annotation = ha
)

pdf(paste0(odir,"GSEA/fgsea_byday_bycelltype.pdf"), width = 18, height = 10)
draw (oh, column_title = "GSEA Old vs Young (REACTOME)",
      heatmap_legend_side = "bottom", padding = unit(c(2, 2, 2, 60), "mm"))
dev.off()

# ========================= HALLMARK : not done ================================
# no needed as not performed

# ===================== Test some other plot =============================
m_ = readRDS(paste0(odir,"GSEA/fgsea_matrices4heatmaps.rds" ))
library(pheatmap)
library(cowplot)
library(heatmap3)
library(RColorBrewer)
library(gridExtra)

names(m_)
names(m_[["D7"]])
# test on single mapcolor with up and down pathways
k = "D4"
CT = "M2"
pdf(paste0(odir,"t.pdf"), width=8, height=5)
heatmap3(as.matrix(m_[[k]][[CT]][["UP"]]),  
         scale = "none",
         col=colorRampPalette(c("gray","firebrick3"))(256), 
         Colv = NA, Rowv = NA,
         cexRow = 0.7,
         cexCol = 0.7,
         ColSideWidth = ncol(m_[[k]][[CT]][["UP"]]),
         margins = c(10,10),
         main = paste(k, CT, "UP"))

heatmap3(as.matrix(m_[[k]][[CT]][["DOWN"]]),  
         scale = "none",
         col=colorRampPalette(c("navy", "gray"))(256), 
         Colv = NA, Rowv = NA,
         cexRow = 0.7,
         cexCol = 0.7,
         ColSideWidth = ncol(m_[[k]][[CT]][["DOWN"]]),
         margins = c(10,10),
         main = paste(k, CT, "DOWN"))
dev.off()
print("end")



