cd ~/BulkAnalysis_plusNetwork/natmiOut_TPM

python ~/NATMI/DiffEdges.py --refFolder OldD0 --targetFolder YoungD0  --interDB lrc2p --weightType mean --out Diff_D0
python ~/NATMI/DiffEdges.py --refFolder OldD2 --targetFolder YoungD2  --interDB lrc2p --weightType mean --out Diff_D2
python ~/NATMI/DiffEdges.py --refFolder OldD4 --targetFolder YoungD4 --interDB lrc2p --weightType mean --out Diff_D4
python ~/NATMI/DiffEdges.py --refFolder OldD7 --targetFolder YoungD7 --interDB lrc2p --weightType mean --out Diff_D7

cd ~/BulkAnalysis_plusNetwork/natmiOut_CountNormalised

python ~/NATMI/DiffEdges.py --refFolder OldD0 --targetFolder YoungD0  --interDB lrc2p --weightType mean --out Diff_D0
python ~/NATMI/DiffEdges.py --refFolder OldD2 --targetFolder YoungD2  --interDB lrc2p --weightType mean --out Diff_D2
python ~/NATMI/DiffEdges.py --refFolder OldD4 --targetFolder YoungD4 --interDB lrc2p --weightType mean --out Diff_D4
python ~/NATMI/DiffEdges.py --refFolder OldD7 --targetFolder YoungD7 --interDB lrc2p --weightType mean --out Diff_D7
