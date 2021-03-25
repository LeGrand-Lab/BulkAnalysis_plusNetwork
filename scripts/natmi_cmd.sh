
cd ~/NATMI

python ExtractEdges.py \
  --species mouse \
  --emFile ~/bulk_analysis/data/meanTPMYoungD7.txt \
  --interDB lrc2p \
  --idType ensembl \
  --out ~/bulk_analysis/natmiOut/Young_D7/ \
  --coreNum 4
