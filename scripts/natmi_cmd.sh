
cd ~/NATMI

for i in $(ls ~/BulkAnalysis_plusNetwork/inDataNatmi/TPM_*); do 
	echo $i | awk -F '[._]' '{print $3}'
	CW=$(echo $i | awk -F '[._]' '{print $3}')
	python ExtractEdges.py \
		--species mouse \
		--emFile ~/BulkAnalysis_plusNetwork/inDataNatmi/TPM_${CW}.txt \
		--annFile ~/BulkAnalysis_plusNetwork/inDataNatmi/annot_${CW}.txt\
		--interDB lrc2p \
		--idType ensembl \
		--out ~/BulkAnalysis_plusNetwork/natmiOut_TPM/${CW}/ \
		--coreNum 8
	echo ""
done

for i in $(ls ~/BulkAnalysis_plusNetwork/inDataNatmi/CountNormalised_*); do 
	echo $i | awk -F '[._]' '{print $3}'
	CW=$(echo $i | awk -F '[._]' '{print $3}')
	python ExtractEdges.py \
		--species mouse \
		--emFile ~/BulkAnalysis_plusNetwork/inDataNatmi/CountNormalised_${CW}.txt \
		--annFile ~/BulkAnalysis_plusNetwork/inDataNatmi/annot_${CW}.txt\
		--interDB lrc2p \
		--idType ensembl \
		--out ~/BulkAnalysis_plusNetwork/natmiOut_CountNormalised/${CW}/ \
		--coreNum 8
	echo ""
done
