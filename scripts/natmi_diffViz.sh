cd ~/BulkAnalysis_plusNetwork/natmiOut/

for i in Diff_*/; do
	python ~/NATMI/VisInteractions.py \
	   --sourceFolder $i \
	   --interDB lrc2p \
	   --detectionThreshold 0.3 \
	   --drawNetwork y \
	   --drawClusterPair y \
	   --plotWidth 12 \
	   --plotHeight 10 \
	   --layout kk \
	   --fontSize 8 \
	   --edgeWidth 0 \
	   --clusterDistance 1
done

#python  --sourceFolder Diff_D7 --inderDB lrc2p --drawClusterPair y

