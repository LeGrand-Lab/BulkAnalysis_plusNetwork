cd ~/BulkAnalysis_plusNetwork/natmiOut

python ~/NATMI/DiffEdges.py --refFolder YoungD0 --targetFolder OldD0 --interDB lrc2p --weightType mean --out Diff_D0
python ~/NATMI/DiffEdges.py --refFolder YoungD2 --targetFolder OldD2 --interDB lrc2p --weightType mean --out Diff_D2
python ~/NATMI/DiffEdges.py --refFolder YoungD4 --targetFolder OldD4 --interDB lrc2p --weightType mean --out Diff_D4
python ~/NATMI/DiffEdges.py --refFolder YoungD7 --targetFolder OldD7 --interDB lrc2p --weightType mean --out Diff_D7

