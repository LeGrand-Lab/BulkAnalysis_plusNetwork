cd ~/Softwares_bio/NATMI

for i in $(ls ~/BulkAnalysis_plusNetwork/NatmiData/inDataNatmi/TPM_*); do 
	echo $i | awk -F '[_.]' '{print $3}'
	CW=$(echo $i | awk -F '[_.]' '{print $3}')
	python VisInteractions.py \
	  --sourceFolder ~/BulkAnalysis_plusNetwork/NatmiData/natmiOut_TPM/${CW}/ \
	  --interDB lrc2p \
	  --weightType mean \
	  --detectionThreshold 0.6 \
	  --plotFormat pdf \
	  --drawNetwork y \
	  --plotWidth 12 \
	  --plotHeight 10 \
	  --layout kk \
	  --fontSize 8 \
	  --edgeWidth 0 --maxClusterSize 0 --clusterDistance 1 
	echo ""
done



for i in $(ls ~/BulkAnalysis_plusNetwork/data/CountNormalised/CountNormalised_*); do 
	echo $i | awk -F '[_.]' '{print $3}'
	CW=$(echo $i | awk -F '[_.]' '{print $3}')
	python VisInteractions.py \
	  --sourceFolder ~/BulkAnalysis_plusNetwork/NatmiData/natmiOut_CountNormalised/${CW}/ \
	  --interDB lrc2p \
	  --weightType mean \
	  --detectionThreshold 0.6 \
	  --plotFormat pdf \
	  --drawNetwork y \
	  --plotWidth 12 \
	  --plotHeight 10 \
	  --layout kk \
	  --fontSize 8 \
	  --edgeWidth 0 --maxClusterSize 0 --clusterDistance 1 
	echo ""
done