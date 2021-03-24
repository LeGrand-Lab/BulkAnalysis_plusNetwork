### what scripts do:
##### Finished:

- `GET_ELEMENTS` : take from count and Stringtie files all information,
yields into 'data/' TPM FPKM Covariance (type csv), and a rds object containing raw counts matrix.

- `get_biotypes_bulk` : connects to biomart, classifies all non-null count genes, yields "data/biotype_bulk.csv".

- `protcod_TPMandCounts` : connects matrices to biotype csv file  to obtain matrices containing only protein coding genes, and also figures into single pdf: 
		- `plotsPrelim/pre_filtering_figures.pdf`
		- `data/protcod_counts.rds`
		- `data/protcod_TPM.rds`

- `matrices_prefiltering` : accepts in input protein-coding counts and TPM, yields outs:
		- "plotsPrelim/" : several PCA
		- "data/prefiltered_counts.rds"
		- "data/prefiltered_TPM.rds"
		- "data/metadata.rds"

##### Unfinished:
- `signaturestypes_Y.R` : 

NOTE:

- `fromhome.R` : 

- `global_stats` : 
- `dyn_oldVSyoung.R` : dynamics **by celltype subset** (Old vs Young across all days)

####### TODO:
find script that creates "gene_diff.pdf"

