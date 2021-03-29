
cd ~/NATMI

for i in $(ls ~/BulkAnalysis_plusNetwork/data/meanTPM_*); do 
	echo $i | awk -F '[._]' '{print $3}'
	CW=$(echo $i | awk -F '[._]' '{print $3}')
	python ExtractEdges.py \
		--species mouse \
		--emFile ~/BulkAnalysis_plusNetwork/data/meanTPM_${CW}.txt \
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
