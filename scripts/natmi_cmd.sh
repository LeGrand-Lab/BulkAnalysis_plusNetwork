
cd ~/Softwares_bio/NATMI

for i in $(ls ~/BulkAnalysis_plusNetwork/NatmiData/inDataNatmi/TPM_*); do 
	echo $i | awk -F '[._]' '{print $3}'
	CW=$(echo $i | awk -F '[._]' '{print $3}')
	python ExtractEdges.py \
		--species mouse \
		--emFile ~/BulkAnalysis_plusNetwork/NatmiData/inDataNatmi/TPM_${CW}.txt \
		--annFile ~/BulkAnalysis_plusNetwork/NatmiData/inDataNatmi/annot_${CW}.txt\
		--interDB lrc2p \
		--idType ensembl \
		--out ~/BulkAnalysis_plusNetwork/NatmiData/natmiOut_TPM/${CW}/ \
		--coreNum 8
	echo ""
done

for i in $(ls ~/BulkAnalysis_plusNetwork/data/CountNormalised/CountNormalised_*); do 
	echo $i | awk -F '[._]' '{print $3}'
	CW=$(echo $i | awk -F '[._]' '{print $3}')
	python ExtractEdges.py \
		--species mouse \
		--emFile ~/BulkAnalysis_plusNetwork/data/CountNormalised/CountNormalised_${CW}.txt \
		--annFile ~/BulkAnalysis_plusNetwork/NatmiData/inDataNatmi/annot_${CW}.txt\
		--interDB lrc2p \
		--idType ensembl \
		--out ~/BulkAnalysis_plusNetwork/NatmiData/natmiOut_CountNormalised/${CW}/ \
		--coreNum 8
	echo ""
done
