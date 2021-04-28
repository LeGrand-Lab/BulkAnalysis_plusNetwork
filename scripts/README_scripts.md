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
	
- `tpm_prep4Tau` : saves into `data/` the 'meanTPM(..).txt' matrices by age 
and day, suitable for specificity index calculation (Tau) .

- `tpm_prep4Nat` : saves into `inDataNatmi/` the TPM_(..).txt and annotation
.txt files suitable for Natmi Ligand-Receptor network generation.

##  Current folder:

- `calc_Tau_Specificity.R` and `calc_Tau_figures.R` are related, the first takes 
TPM matrices for Tau tissue specificity index calculation. The second performs 
related plots, all outputs saved into `Tau/` folder.

- `dynamics_intra_Spec.R` : 
ONLY FOR TISSUE SPECIFIC GENES, analyzes gene expression dynamics  existing 'intra condition', i.e. separately
for 'young' and 'old' age groups. Uses pairwise DESeq2 approach. Results saved into folder `dynamicsIntra_Spec/`.

- `dynamics_intra_doGO.R` :  Takes results from previous script, and performs enrichment
by means of gprofiler2. Preferred terms are Reactome (REACT) and GO. Results saved into
`dynamicsIntra_Spec/Gprofiler2_res/`.

- `dynamics_intra_edgeR_vrs.R`:  for both specific and housekeeping, i.e.  no previous filtering on counts matrix, analyzes global intra group expression dynamics with edgeR following same contrasts as done previously with DESeq2. `dynamics_intra_edgergoKegg.R` performs GO and KEGG enrichment, plots by `dynamics_intra_edgerkeggplots.R` . 
All results saved into `dynintra_edger_extended/`. Too extensive for being interpretable. 

- `exam_INTER_cond_dyn.R` : analyzes if dynamics in expression is differential between Young and Old conditions.
Results saved into `exam_INTER_conditions/dynamic`.

- `exam_Inter_cond_sta.R` : performs classic DE test (Old vs Young) time point by time point (like a "snapshot" of  expression differences at each time point). Results saved into `exam_INTER_conditions/static`

- natmi_cmd.sh and natmi_edges.sh: ran in this order, yield  L-R networks into
natmiOut/ folder.

- pathView_test.R : a simple test with pathView, results into 'path_classic_tools/'.

- yieldsArbitraryKmeansPlot.R : a demo to show that, for "specific" genes,
increasing 'k' yields repetitive patterns (a justification to use 
silhouette testing up to k=8) after kmeans calc. Saves plots into dynamicsIntra/.


JohaGL 2021


