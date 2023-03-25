## BulkAnalysis_plusNetwork

This repository contains a bio-informatic analysis related to Skeletal Muscle Repair in Aged vs Juvenile backgrounds RNA-seq time-series (Chazaud + Feige + LeGrand teams).
,

RNA-seq data are pre-processed with nf-core pipeline https://nf-co.re/rnaseq (Mapping with STAR, counting with FeatureCount and StringTie)
,

`scripts` folder contains R code which is detailed in [this relative link](scripts/README_scripts.md). 
In summary:
- they run DEG and GSEA analyses with result in exam_INTER_conditions/static/ and reports are in https://github.com/LeGrand-Lab/Ageing-impact_in_gene_expression_on_skeletal_muscle_repair.
- Prepare data for ligand-receptor network and finally launch NATMI  https://github.com/asrhou/NATMI. 

Many thanks to Dr.Rui Hou for advice concerning NATMI to use it specifically in our context.

### Ligands-Receptors part

The dedicated analysis started as shown in `networks_explore/` folder: the "Edges" tables of the Natmi results were processed by the jupyter notebook 1_  to generate the .ml files. Those .ml files were then copied into the Shiny app  `AgeingMice_shinyApp`. This application allow to navigate, by day, the results of DEGs, L-R networks and enrichment. It can be deployed locally, see "Technical note" below. 

####  Technical note

For conda users we provide the two .yml files  [here](scripts/installs/muscle_conda_yml_envs) to create, respectively, the Python and R environments with all the libraries and packages (the part "networks_explore" requires the Python environment, whereas  the Shiny app requires the R environment). 

Fast run (after conda steps) : 
```
cd $HOME/BulkAnalysis_plusNetwork
R -e "shiny::runApp('AgeingMice_shinyApp')"
```




Otherwise, you must install the packages in your preferred way and run the app in RStudio.



###  More about Pathways in this dataset

See also 
https://github.com/LeGrand-Lab/Ageing-impact_in_gene_expression_on_skeletal_muscle_repair (authored by Pauline) for interactive reports regarding advanced reactome pathways analysis on full dataset, and more!


##### Contact info

* juana7@gmail.com (Johanna), 

* pauline.moulle@cnrs.fr

* william.jarassier@univ-lyon1.fr


### authors
[johaGL](https://github.com/johaGL/)
M2 intern Bio-informatics
Université Claude Bernard Lyon 1
Institut Neuromyogène INMG

[PaulineMoulle](https://github.com/PaulineMoulle)
Bioinformatique Engineer
Institut Neuromyogène INMG
