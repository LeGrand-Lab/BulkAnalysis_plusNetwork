# Recap,  methods  and projections for future

This document explain the files at `exam_INTER_conditions/static/`. Scripts are located in `scripts/` as usual.

Why are we interested in "static" approach (and not dynamic one) ?. The answer:

1. **Neutro** cell type is only present at day D2, if we only take account of dynamic approach, we will neglect this important celltype for differences Old vs Young at D2 (furthermore, it has interactions with other celltypes at D2, see numeral 3).
2. We are interested in seing what happens to Pathways across time, with details into every single day, so we need to run Pathway analysis day by day.
3. We need to integrate DE information (FC==log2FoldChange, padj, etc)  to Ligand-Receptor networks. These Networks were built **by day**.

## Recap

### A. DEGs list / DE full information: 

Under the assertion that Young is the "control" scenario, a classical DE test, Old vs Young was performed by cell type by day. Full results (not DEGs but entire metrics for all tested genes) were saved as a 320346 x 9 dataframe, compressed into a rds file `rds/shot_rds_full.rds`. Softfiltered and filtered versions (the DEGs lists) are available at `csv/`, for sharing with biologists. For bioinformaticians, use full rds allways, and set filters (or not) depending on the needs. 

We found:
- some strong effects that were not significant, which can be in part explained by few biological replicates (low statistical power)
- some weak effects that were significant (where "Old" not old enough, so the differences compared to Young are weak ?)
- some strong and significant effects, mostly present in FAPs and sCs. This is somewhat expected by biologists.
- A lot of genes not being differentially expressed among ages.

### B. Tau indexes:

johaGL calculated (from TPM matrices) Tau tissue specificity indexes for each **day+age**. She determined with a simple 'which.max' the most likely celltype for which a gene was a specific marker. 

 The goal is to **integrate this information to full DE metrics** (as explained in section A.). One single Tau value per gene is therefore required.

**Tau classification and consensus** : 
When checking this confusion matrix (by johaGL):

|	         |high:old	|low:old	| absent:old	|
|--------------|--------------|--------------|--------------|
|high:young  | TrueSpecific	|Specific downregulated-old 	|  Specific absent-old |
|low:young   | Housekeeping Upregulated-old	|TrueHousekeeping | Housekeeping absent-old |
|absent:young | Exclusive Old | Exclusive Old | Not Available  |


One can think that taking **Young** Tau values would be enough (they differentiate specific or housekeeping categories). However, the *worst case scenario* is when, for a gene, Tau does not exist (null TPM) in Young, but it only exists in Old,  which means it has exclusive expression in Old (from TPM values). To handle this case, script `calc_Tau_Specificity.R` creates a consensus dataframe with a column "exclusiveOld": For a *i* gene, let *yv_i* be the young value, *ov_i* the old value. If *yv_i* is NA, then add to table the *ov_i* value, and set "exclusiveOld" as 1. Otherwise, set 0 and add *yv_i*.  The *ages consensus*  was built as a list of dataframes (one dataframe by **day** ) :  `conseTau_ensemblid.rds` (`Tau/` dir).

Note:  Even for those genes whom Tau was low, there is a "whichMax" (celltype), so be careful.

Note about **cutoff** for defining "Specific" category: Initially, tau > 0.7 was fixed. However in M1 at D2, Ly96 has tau 0.319711138642448, but Ly96 is macrophage specific. Other example in same M1 (D2) is Cd68, with tau 0.428672631480618. So for fixing a cutoff for future applications, johaGL recommends 0.3.

### C. GSEA  

1. **GSEA by day by cell type**

Mixed housekeeping and specific (Tau info not integrated). Full DE info dataframes were splitted into up and down datasets, then ordered by descending absFC and smallest to biggest padj, cutting at first 500 genes (500up + 500down = 1000 genes).  The resulting "blinded" cutoffs are:

|  day_celltype | inputsize | maxpadj | minabslfc |
|  ------------- | ------------ | ----------- | ---------- |
|   D0_ECs | 1000 | 0.833204755812667 | 0.00431642344880973  |
|   D0_FAPs | 1000 | 0.211068949639519 | 0.0218738426298392  |
|   D0_sCs | 1000 | 0.200496565670897 | 0.0328251937026974  |
|   D2_ECs | 1000 | 0.99907937632229 | 5.86061346118547e-07  |
|   D2_FAPs | 1000 | 0.445153360596376 | 0.00860019026213538  |
|   D2_M1 | 1000 | 0.999797985224188 | 3.27646804795741e-09  |
|   D2_M2 | 1000 | 0.656706752930478 | 0.0138428791183234  |
|   D2_Neutro | 1000 | 0.445119625266764 | 0.0322225565883877  |
|   D2_sCs | 1000 | 0.0141307968995322 | 0.405389381381288  |
|   D4_ECs | 1000 | 0.308373315453155 | 0.0176139001908065  |
|   D4_FAPs | 1000 | 0.393112574721973 | 0.014523666546943  |
|   D4_M1 | 1000 | 0.913181448464194 | 1.07783951528171e-05  |
|   D4_M2 | 1000 | 0.83325821418217 | 0.00177164842691496  |
|   D4_sCs | 1000 | 0.034332816186047 | 0.176309694233418  |
|   D7_ECs | 1000 | 0.18952286163962 | 0.0214819598200951  |
|   D7_FAPs | 1000 | 0.308291979666569 | 0.0140952915537856  |
|   D7_M2 | 1000 | 0.712679100846678 | 0.00874214071918688  |
|   D7_sCs | 1000 | 0.202909455037087 | 0.0284891378559077  |


These cutoffs are the only possible to keep 1000 genes in input, as fgsea padj results improved with input size in our dataset. Moreover, genes showing weak effects (or null) constitute a spectrum that does not deformate results, because strong effects will dominate and produce bigger absolute NES.
The smallest FDR (padj) values obtained across GSEA results are :  

```
            ECs         FAPs           M1          M2    Neutro          sCs
D0 9.974770e-01 1.681250e-09           NA          NA        NA 1.634265e-02
D2 6.139028e-02 2.070819e-05 0.0000001566 0.001139709 0.6745883 1.180000e-08
D4 2.622367e-06 8.362185e-01 0.7402079002 0.055212297        NA 1.133333e-08
D7 1.054715e-02 5.747208e-01           NA 0.294533841        NA 7.095938e-06
```

2. **Try to incorporate Tau to DE results by day by celltype**

This time the aim was to "integrate" Tau information to DE information. By Day by Cell type.

For each produced *mixed* dataframe, there were some discoveries:

- rows containing  NA padj values were sometimes associated with Tau > 0.5, but as expected FC was near zero.
```
mix %>% filter(is.na(padj) & !is.na(Tau) & Tau > 0.5) %>% slice_max(abs(log2FoldChange),n=1)
# A tibble: 1 ?? 15
  baseMean log2FoldChange lfcSE  pvalue  padj id                 day   type  symbol.x symbol.y Tau               class    whichMAX nbMAX exclusiveOld
     <dbl>          <dbl> <dbl>   <dbl> <dbl> <chr>              <chr> <chr> <chr>    <chr>    <chr>             <chr>    <chr>    <chr> <chr>       
1     18.8        -0.0616 0.139 0.00679    NA ENSMUSG00000030559 D7    sCs   Rab38    Rab38    0.979686856413078 specific sCs      1     0      
```

- Thousands of genes appearing in DE info do not have Tau values, potential explanations: 
	* DE comes from count matrices whereas Tau comes from TPM matrices. 
	* Calculating Tau implies to bring all libraries (for a given day) together. DE analysis, on the other side, was performed on separated libraries. In consequence, matching 'whichMAX' (in Tau info) with 'type' (in DE info) does not guarantees true correspondance. 

- Very few of the genes showed any effect (absFC > 1), so there was no utility of filtering by Tau. In conclusion, for contrast performed *Old vs Young by day by cell type*, using Tau index has no sense.

In any case, when discussing with biologists it was concluded that it would be very interesting to filter by Tau values, which give an approximation to cell-type's markers and their behaviour accros time. But in this case, another contrast is necessary (Old vs Young by day).


3. **GSEA by day (gathered celltypes)** 

It is performed by script :
  `exam_Intx_TauSelec.R` .
  
See 'exam_Intx_report.pdf' file to see if matrices filtered by Tau are suitable for DESeq2 (the answer is yes, dispersions are coherent with the model required by DESeq2).

## Projection:
To make integrated heatmaps, or network "graph" style representations integrating the notion of ages differences across time, and if possible, add Reactome pathways terms in some way. Hard work yet to come ! 

## Useful sources:
- [https://bioinformatics-core-shared-training.github.io/RNAseq_May_2020_remote/html/06_Gene_set_testing.html#fgsea](https://bioinformatics-core-shared-training.github.io/RNAseq_May_2020_remote/html/06_Gene_set_testing.html#fgsea)

 
