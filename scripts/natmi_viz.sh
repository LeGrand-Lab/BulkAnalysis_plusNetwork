
cd ~/NATMI

for i in $(ls ~/BulkAnalysis_plusNetwork/data/meanTPM_*); do 
	echo $i | awk -F '[_.]' '{print $3}'
	CW=$(echo $i | awk -F '[_.]' '{print $3}')
	python VisInteractions.py \
	  --sourceFolder ~/BulkAnalysis_plusNetwork/natmiOut/${CW}/ \
	  --interDB lrc2p \
	  --weightType mean \
	  --detectionThreshold 0.2 \
	  --plotFormat pdf \
	  --drawNetwork y \
	  --plotWidth 12 \
	  --plotHeight 10 \
	  --layout kk \
	  --fontSize 8 \
	  --edgeWidth 0 --maxClusterSize 0 --clusterDistance 1 
	echo ""
done



#python VisInteractions.py \
#  --sourceFolder ~/bulk_analysis/natmiOut/Young_D7/ \
# --interDB lrc2p \
# --weightType mean \
#  --detectionThreshold 0.2 \
#  --plotFormat pdf \
#  --drawNetwork y \
#  --plotWidth 12 \
#  --plotHeight 10 \
#  --layout kk \
#  --fontSize 8 \
#  --edgeWidth 0 --maxClusterSize 0 --clusterDistance 1 
