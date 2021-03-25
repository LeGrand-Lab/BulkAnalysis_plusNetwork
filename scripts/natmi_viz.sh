cd ~/NATMI

python VisInteractions.py \
  --sourceFolder ~/bulk_analysis/natmiOut/Young_D7/ \
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
