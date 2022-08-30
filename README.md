## BulkAnalysis_plusNetwork

This repository contains a bio-informatic analysis related to Skeletal Muscle Repair in Aged vs Juvenile backgrounds RNA-seq time-series (Chazaud + Feige + LeGrand teams).
,

RNA-seq data are pre-processed with nf-core pipeline https://nf-co.re/rnaseq (Mapping with STAR, counting with FeatureCount and StringTie)
,

`scripts` folder contains R code which is detailed in [this relative link](scripts/README_scripts.md). 
In summary:
- they run DEG and GSEA analyses with result in exam_INTER_conditions/static/ and reports in reports/.
- Prepare data for ligand-receptor network and finally launch NATMI  https://github.com/asrhou/NATMI. 

Many thanks to Dr.Rui Hou for advice concerning NATMI to use it specifically in our context.

Three jupyter notebooks are available in folder `networks_explore`:
* LR_networks_all.ipynb : networks data manipulation, LRinfos and Networkx objects saved into dictionnaries, then into two binary files.
* LR_nx_to_igraph.ipynb : transforms Networkx objects into igraph objects, then into separate graphml files
* Young_D7_test_initial.ipynb : checks properties and how to get connected components (using Young D7 for demo purposes)

Folder `AgeingMice_shinyApp` contains Shiny application code, in active development.


For more information write to:
deisy.lascroux@etu.univ-lyon1.fr (or juana7@gmail.com), 
william.jarassier@univ-lyon1.fr
pauline.moulle@cnrs.fr

### author
[johaGL](https://github.com/johaGL/)
M2 intern Bio-informatics
Université Claude Bernard Lyon 1
Institut Neuromyogène INMG

[PaulineMoulle](https://github.com/PaulineMoulle)
Bioinformatique Engineer
Institut Neuromyogène INMG
