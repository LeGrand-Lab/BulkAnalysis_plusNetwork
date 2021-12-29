README

Reference dataset: /home/bioinfo/BulkAnalysis_plusNetwork/natmiOut/YoungD0
Target dataset: /home/bioinfo/BulkAnalysis_plusNetwork/natmiOut/OldD0

Reference dataset has 10382 edges
Target dataset has 10027 edges

cluster_comparison.xlsx: cluster size and fraction of each cluster in two datasets
cluster_size_comparison.pdf: comparison of the cluster sizes of each cluster in two datasets
cluster_fraction_comparison.pdf: comparison of the cluster fractions of each cluster in two datasets

Delta_edges_xxx folder: variations in the edges of two datasets using the certain ligand-receptor pair list
	|->Disappeared_xxx.csv: edges that only detected in the reference dataset
	|->Appeared_xxx.csv: edges that only detected in the target dataset
	|->UP-regulated_xxx.csv: edges that have higher expression weights in the target dataset
	|->DOWN-regulated_xxx.csv: edges that have lower expression weights in the target dataset
	|->Stable_xxx.csv: edges that have identical expression weights in two datasets
	|->All_edges_xxx.csv: all edges that detected in both dataset

Sending cluster: clusters that expresses the ligand
Target cluster: clusters that expresses the receptor
Ligand/receptor symbol: official gene symbol of the detected ligand/receptor
Delta ligand/receptor expression: change of ligand/receptor expression values in two datasets
Delta ligand/receptor specificity: change of ligand/receptor expression value derived specificities in two datasets
Delta edge expression weight: change of the edge expression weights in two datasets
Delta edge specificity weight: change of the edge specificity weights in two datasets
Log2-transformed fold change of edge expression weight: ratio of the edge expression weight in the target dataset to that in the reference dataset
Log2-transformed fold change of edge specificity weight: ratio of the edge specificity weight in the target dataset to that in the reference dataset

