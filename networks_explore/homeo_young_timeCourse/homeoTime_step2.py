#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
- extract top variable ligands (tunable n)
- find all receptors linked to top variable ligands from step 1 (L-R pairs)
- for each L-R pair, extract a list of edge weights accross days (4 values):
     use column 'Edge average expression weight' by each day;
     group by ligand
- for each ligand, calculate variance of edge weights across time
- pick top variable edges (tunable n)                
- plot and write csv

Created on Wed Oct 13 11:54:37 2021

@author: johaGL
"""

import os
import pickle
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from scipy import stats
import copy
import seaborn as sns
from boxFunClass import *

def filldaywvectorstodict(topligdfrec, dicow, day):
    """
    output: { ('Adam10_Neutro', 'Ehpa3_sCs') : [w1,w2,w3,w4] ...}
    """
    weis_ = [-1,-1,-1,-1]    
    dki = {'D0' : 0, 'D2' : 1, 'D4' : 2, 'D7' : 3}  # indexes in vector 
    for i, row in topligdfrec.iterrows():
        fromvertex = row['uniq_Ligand_symbol']
        tovertex = row['uniq_Receptor_symbol']
        weiday = row['Edge average expression weight']
        try:
            dicow[(fromvertex,tovertex)][dki[day]] = float(weiday)
        except KeyError:
            dicow[(fromvertex,tovertex)] = copy.deepcopy(weis_)
            dicow[(fromvertex,tovertex)][dki[day]] = float(weiday)  
    return dicow

    
def reorderdico_fillvariancew(dicow):
    """
    output: dico { ligand : [(rec1, varianceweight, [w1,w2,w3,w4]),(rec2, .. }
    Note:  tup[0] is receptor, tup[1] is varianceweight, tup[2] the vector
    """
    varsdico = {}
    for tupk in dicow.keys():        
        ligand, receptor = tupk[0], tupk[1]
        truewei = [i for i in dicow[tupk] if i > -1]
        if len(truewei) > 1:
            v = np.var(truewei)
        elif len(truewei) == 1:
            v = truewei[0]
        try :
            varsdico[ligand].append((receptor, v, dicow[tupk] ))
        except KeyError:
            varsdico[ligand] = [(receptor, v, dicow[tupk])]
    return varsdico

########## prepare variables
lr = pickle.load( open( '../graphobjs/dictio_lr.p', 'rb') )
hom = lr['Young']
days = ['D0', 'D2', 'D4', 'D7']
allcelltypes =  [ "ECs" ,
    "FAPs"  ,
    "Inflammatory-Mac"  ,
    "Resolving-Mac"  ,
    "Neutrophils" ,
    "MuSCs" ]

#fileligs = "selectedLigands.txt"
fileligs = "selectedLigandsVersion3.txt"
celltycolors = importcelltycolorsdico()
Nligs = 15
NtoprecBylig = 4
print(f"'\n'Finding top {Nligs} ligands, every celltype, you can modify this number")
topdicoligs = yieldtopUnique(fileligs, allcelltypes, Nligs)
print("printing top Satellite Cells to show up what we got")
print(topdicoligs['MuSCs'])
topligslt_ = []
for k in topdicoligs.keys():
    for tup in topdicoligs[k]: 
        topligslt_.append('_'.join([tup[0],k])) # ex: 'Ptn_FAPs'
        

######### find only for muscular satellites the receptors for selected ligands
print("Show ligands from all cell types targeting receptor in Satellite Cells")
prepsatelly = {}
for day in days:
    daydf = hom[day].frame
    prepadf = daydf[daydf['Target cluster'] == 'MuSCs']
    topligdfrec = prepadf[prepadf['uniq_Ligand_symbol'].isin(topligslt_)]
    prepsatelly = filldaywvectorstodict(topligdfrec, prepsatelly,day)
    
defisatellite = sortandfilterltupdic(reorderdico_fillvariancew(prepsatelly), 1, True)  
satellireceptors = set()
liggs = set()
    
for k in defisatellite.keys():
    lig = k
    liggs.add(k)
    listofreceptors = defisatellite[k]
    for tup in listofreceptors:        
        rec = tup[0]
        satellireceptors.add(rec)
    
# order ligands by celltype origin, alphabetical
preligorder = Sort([i.split("_") for i in liggs], 1, False)
ligorder = ["_".join(l_) for l_ in preligorder] 

"""
Heatmap with weights across time
"""
#print(ligorder)
#print(satellireceptors)
rowlabels = [] 
symboligs = [] ; symborecs = []
origcelltype = []
targetcelltype = [] 
day0 = [] ; day2 = [] ; day4 = [] ; day7 = []

for k in ligorder:
    tmpli = k.split("_")
    lig, origct = tmpli[0], tmpli[1]
    listofreceptors = defisatellite[k] # use key directly to catch lig infos
    for tup in listofreceptors:  # tup
        tmpre = tup[0]
        tmpre = tmpre.split("_")
        rec, targrec = tmpre[0], tmpre[1]
        rowlabels.append(" --> ".join([lig,rec]))
        symboligs.append(lig)
        symborecs.append(rec)
        origcelltype.append(origct)
        targetcelltype.append(targrec)
        day0.append(tup[2][0])
        day2.append(tup[2][1])
        day4.append(tup[2][2])
        day7.append(tup[2][3])
        
dfplot = pd.DataFrame({  'LR_pairs': rowlabels,
                       '0' : np.log10(day0),
                        '2' : np.log10(day2),
                       '4' : np.log10(day4),
                       '7' : np.log10(day7)})
dfplot = dfplot.set_index(['LR_pairs'])


print(f'\nprinting celltypes colors dictionary (celltycolors)')
print(celltycolors)
row_co = pd.DataFrame({"Sender_Type" : [celltycolors[ct] for ct in origcelltype],
                      "LR_pairs" : rowlabels})
row_co = row_co.set_index(['LR_pairs'])
g = sns.clustermap(dfplot, center=1,
                    row_cluster=False, 
                    col_cluster=False,
                    row_colors=row_co,
                    cmap="YlGnBu", 
                    yticklabels=True, 
                    xticklabels=True,
                    cbar_pos=(0.2, 0.83 ,  0.025, 0.08),
                    cbar_kws={'orientation': 'vertical',
                              'label':'log10 edge weight'  },
                    figsize=(8,18)
                    )

for label in celltycolors.keys():
    print(label)
    g.ax_col_dendrogram.bar(1.5, 0 , color=celltycolors[label],
                            label=label, linewidth=0)

g.ax_col_dendrogram.legend(loc="center right", ncol=2, title="Sender cell types")
ax = g.ax_heatmap
ax.set_xlabel("DAYS")
ax.tick_params(bottom=False) 
#g.cax.set_position([.03, .2, .03, .45])
g.fig.suptitle("Ligand-Receptor pairs: \n top signals to Satellite Cell receptors" ,
                    fontsize=24, horizontalalignment='center' )
g.savefig("MainSignalstoSatelliteCells.pdf", pad_inches=1)
plt.show()

#### save a csv file with data used to plot 'MainSignalstoSatelliteCells.pdf'

satellidf = pd.DataFrame({'Sender group' : origcelltype,
                       'Ligand' : symboligs,
                       'Receptor' : symborecs,
                       'Target group' : targetcelltype ,
                       'weight D0' : day0,
                        'weight D2' : day2,
                        'weight D4' : day4,
                        'weight D7' : day7})

satellidf.to_csv("MainSignalstoSatelliteCells.csv", sep=",")


############## find for all celltypes the receptors for selected ligands        
print(f"'\n'Getting the receptors for selected topligs list of tuples (topligslt_)") 
prepall = {} 
for day in days:
    daydf = hom[day].frame 
    topligdfrec = daydf[daydf['uniq_Ligand_symbol'].isin(topligslt_)]
    #print(topligdfrec[['uniq_Receptor_symbol', 'Edge average expression weight', 'uniq_Ligand_symbol']])
    prepall = filldaywvectorstodict(topligdfrec, prepall, day)         

print(f"\n Selecting top {NtoprecBylig} receptors for each of the selected ligands")
defdicofiltered = sortandfilterltupdic(reorderdico_fillvariancew(prepall), NtoprecBylig, True)

############## APENDIX

def doasmatrix(defisatellite, ligorder, satellireceptors, day=None):
    """
     Matrix plot idea :  with signals that arrive to satellite 
     I saw not a good idea: very sparse, time not easy to show in same plot
    """     
    def createdicoindexes(alistofstrings):
        # map of indexes, matches with matrix indexes
        d = {}
        inde = 0
        while inde < len(alistofstrings):
            d[alistofstrings[inde]] = inde
            inde += 1
        return d
    if day is None:
        day = 2
    ligosin = createdicoindexes(ligorder)
    recosin = createdicoindexes(list(satellireceptors))
    matrice = np.empty(shape=(len(ligorder), len(satellireceptors)))
    for k in defisatellite.keys():
        lig = k
        i = ligosin[lig]
        listofreceptors = defisatellite[k]
        for tup in listofreceptors:  # tup
            rec = tup[0]
            j = recosin[rec]
            matrice[i,j] = tup[2][0] #tup[2]  is the 4 weights vector    
    
    ax = sns.heatmap(matrice, cmap="Greens")
    plt.show()
    return 0
doasmatrix(defisatellite, ligorder, satellireceptors, day=2)
print("BAD, very sparse matrix, not the good choice; Time not plottable in 1 step")






  


