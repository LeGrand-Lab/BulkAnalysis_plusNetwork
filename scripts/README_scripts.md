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
 it does not filter any gene, kepts good quality libraries and yields outs:
	* "data/prefiltered_counts.rds"
	* "data/prefiltered_TPM.rds"
	* "data/metadata.rds"
	* "plotsPrelim/" : several PCA, and a complex pheatmap
.	
	
- `tpm_prep4TauNat` : saves into `data/` the .txt matrices by age and day, 
suitable for specificity index calculation (Tau) and Natmi Ligand Receptor 
nework creation.

##  Current folder:

- `calc_Tau_Specificity.R` and `calc_Tau_figures.R` are related, the first takes 
TPM matrices for Tau tissue specificity index calculation. The second performs 
related plots, all outputs saved into `Tau/` folder

- `dynamics_intra.R` : 
analyzes gene expression dynamics existing 'intra condition', i.e. separately
for 'young' and 'old' age groups. Results saved into folder `dynamicsIntra/`.

- `


- `;;;.R` : 

- `;;;` : 
- `dyn_oldVSyoung.R` : dynamics **by celltype subset** (Old vs Young across all days)


JohaGL 2021


