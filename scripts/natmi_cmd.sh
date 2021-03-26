
cd ~/NATMI

for i in ~/BulkAnalysis_plusNetwork/data/meanTPM-*; do 
	echo $i | awk -F '[-.]' '{print $2}'
	CW=$(echo $i | awk -F '[-.]' '{print $2}')
	python ExtractEdges.py \
		--species mouse \
		--emFile ~/BulkAnalysis_plusNetwork/data/meanTPM-${CW}.txt \
		--interDB lrc2p \
		--idType ensembl \
		--out ~/BulkAnalysis_plusNetwork/natmiOut/${CW}/ \
		--coreNum 8
	echo ""
done

#python ExtractEdges.py \
#  --species mouse \
#  --emFile ~/bulk_analysis/data/meanTPMYoungD7.txt \
#  --interDB lrc2p \
#  --idType ensembl \
#  --out ~/bulk_analysis/natmiOut/Young_D7/ \
#  --coreNum 4
