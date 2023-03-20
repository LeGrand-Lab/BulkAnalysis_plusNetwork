cd ~/BulkAnalysis_plusNetwork/NatmiData/natmiOut_TPM

python ~/Softwares_bio/NATMI/DiffEdges.py --refFolder OldD0 --targetFolder YoungD0  --interDB lrc2p --weightType mean --out Diff_D0
python ~/Softwares_bio/NATMI/DiffEdges.py --refFolder OldD2 --targetFolder YoungD2  --interDB lrc2p --weightType mean --out Diff_D2
python ~/Softwares_bio/NATMI/DiffEdges.py --refFolder OldD4 --targetFolder YoungD4 --interDB lrc2p --weightType mean --out Diff_D4
python ~/Softwares_bio/NATMI/DiffEdges.py --refFolder OldD7 --targetFolder YoungD7 --interDB lrc2p --weightType mean --out Diff_D7

cd ~/BulkAnalysis_plusNetwork/NatmiData/natmiOut_CountNormalised

python ~/Softwares_bio/NATMI/DiffEdges.py --refFolder OldD0 --targetFolder YoungD0  --interDB lrc2p --weightType mean --out Diff_D0
python ~/Softwares_bio/NATMI/DiffEdges.py --refFolder OldD2 --targetFolder YoungD2  --interDB lrc2p --weightType mean --out Diff_D2
python ~/Softwares_bio/NATMI/DiffEdges.py --refFolder OldD4 --targetFolder YoungD4 --interDB lrc2p --weightType mean --out Diff_D4
python ~/Softwares_bio/NATMI/DiffEdges.py --refFolder OldD7 --targetFolder YoungD7 --interDB lrc2p --weightType mean --out Diff_D7
