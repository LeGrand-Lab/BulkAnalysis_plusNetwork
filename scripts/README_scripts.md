## what these local scripts do :
###  `prepareData` folder:

They are all R scripts, in order of execution:

- `get_expr_to_csv_rds` : take from count and Stringtie files all information,
yields into 'data/' TPM FPKM Covariance (type csv), batchesinfo.csv, 
and a rds object containing raw counts matrix.

- `get_biotypes_bulk` : connects to biomart, classifies all non-null count genes, 
yields "data/biotype_bulk.csv".

- `get_protcod_TPMandCounts` : connects matrices to biotype csv file  
to obtain matrices containing only protein coding genes, and also figures into 
pdf files : 
	* `data/protcod_counts.rds`
	* `data/protcod_TPM.rds`
	* `plotsPrelim/`: biotypes.pdf, rawcounts_obs.pdf and TPM_obs.pdf
.

- `libraries_kept_andPCA` : accepts in input protein-coding counts and TPM,
 drops out bad quality libraries as well as 'all zero count rows'. 
 Yields *'prefiltered'* outs:
	* "data/prefiltered_counts.rds"
	* "data/prefiltered_TPM.rds"
	* "data/metadata.rds"
	* "plotsPrelim/": several PCA (on vst counts),Spearman's correlation pheatmap,
	.
.	
	
- `tpm_prep4TauNat` : saves into `data/` the 'meanTPM(..).txt' matrices by age 
and day, suitable for specificity index calculation (Tau) and Natmi Ligand Receptor nework creation.

##  Current folder:

- `calc_Tau_Specificity.R` and `calc_Tau_figures.R` are related, the first takes 
TPM matrices for Tau tissue specificity index calculation. The second performs 
related plots, all outputs saved into `Tau/` folder

- `dynamics_intra.R` : 
analyzes gene expression dynamics existing 'intra condition', i.e. separately
for 'young' and 'old' age groups. Uses pairwise DESeq2 approach. Results saved into folder `dynamicsIntra/`.

- `


- `;;;.R` : 

dynamics_intra_edgergoKegg.R
dynamics_intra_edgeR_vrs.R
dynamics_intra.R
dyn_inter_oldVSyoung.R
natmi_cmd.sh
natmi_treat_results.R
natmi_viz.sh
path_test.R ===> ??
plotsTimeWiseDyn.R
plotsTimewiseDynVB.R
snapshots_exprs.R


JohaGL 2021


