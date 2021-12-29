## BulkAnalysis_plusNetwork

This repository contains a bio-informatic analysis related to Skeletal Muscle Repair in Aged vs Juvenile backgrounds RNA-seq time-series (Chazaud + Feige + LeGrand teams).

This work is still under development.

`scripts` folder contains R code which is detailed in [this relative link](scripts/README_scripts.md). In summary, they run DE analyses, Tau specificity index on whole protein coding elements, prepare data for ligand-receptor network and finally launch NATMI  https://github.com/asrhou/NATMI. 

Many thanks to Dr.Rui Hou for advice concerning NATMI to use it specifically in our context.

Three jupyter notebooks are available in folder `networks_explore`:
* LR_networks_all.ipynb : networks data manipulation, LRinfos and Networkx objects saved into dictionnaries, then into two binary files.
* LR_nx_to_igraph.ipynb : transforms Networkx objects into igraph objects, then into separate graphml files
* Young_D7_test_initial.ipynb : checks properties and how to get connected components (using Young D7 for demo purposes)

Folder `networks_explore/miniapp_3` contains Shiny application code, in active development.

We will add a Creative Commons licence after submission of our manuscript (very soon).

For more information write to:
deisy.lascroux@etu.univ-lyon1.fr (or juana7@gmail.com), 
william.jarassier@univ-lyon1.fr

### author
[johaGL](https://github.com/johaGL/)
M2 intern Bio-informatics
Université Claude Bernard Lyon 1
Institut Neuromyogène INMG

[PaulineMoulle](https://github.com/PaulineMoulle)
Bioinformatique Engineer
Institut Neuromyogène INMG