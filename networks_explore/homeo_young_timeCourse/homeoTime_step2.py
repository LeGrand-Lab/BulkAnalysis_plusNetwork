#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Oct 13 11:54:37 2021

@author: bioinfo
"""

import os
import pickle
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from scipy import stats
import copy
import seaborn as sns
from boxFunClassVersion3 import *

def filldaywvectorstodict(topligdfrec, deco, day):
    """
    output: { ('Adam10_Neutro', 'Ehpa3_sCs') : [w1,w2,w3,w4] ...}
    """
    weis_ = [-1,-1,-1,-1]    
    dki = {'D0':0,'D2':1, 'D4':2, 'D7':3}  
    for i, row in topligdfrec.iterrows():
        fromvertex = row['uniq_Ligand_symbol']
        tovertex = row['uniq_Receptor_symbol']
        weiday = row['Edge average expression weight']
        try:
            deco[(fromvertex,tovertex)][dki[day]] = float(weiday)
        except KeyError:
            deco[(fromvertex,tovertex)] = copy.deepcopy(weis_)
            deco[(fromvertex,tovertex)][dki[day]] = float(weiday)  
    return deco

    
def reorderdico_fillvariancew(deco):
    """
    output: dico { ligand : [(rec1, varianceweight, [w1,w2,w3,w4]),(rec2, .. }
    Note:  tup[0] is receptor, tup[1] is varianceweight, tup[2] the vector
    """
    varsdico = {}
    for tupk in deco.keys():        
        ligand, receptor = tupk[0], tupk[1]
        truewei = [i for i in deco[tupk] if i > -1]
        if len(truewei) > 1:
            v = np.var(truewei)
        elif len(truewei) == 1:
            v = truewei[0]
        try :
            varsdico[ligand].append((receptor, v, deco[tupk] ))
        except KeyError:
            varsdico[ligand] = [(receptor, v, deco[tupk])]
    return varsdico

############## find for all celltypes the receptors for selected ligands

lr = pickle.load( open( '../graphobjs/dictio_lr.p', 'rb') )
hom = lr['Young']
days = ['D0', 'D2', 'D4', 'D7']
allcelltypes = ['ECs', 'FAPs','M1','M2', 'Neutro', 'sCs']

#fileligs = "selectedLigands.txt"
fileligs = "selectedLigandsVersion3.txt"

Nligs = 15
NtoprecBylig = 4
print(f"'\n'Finding top {Nligs} ligands, every celltype, you can modify this number")
topdicoligs = yieldtopUnique(fileligs, allcelltypes, Nligs)
print("printing top sCs to show up what we got")
print(topdicoligs['sCs'])
topligslt_ = []
for k in topdicoligs.keys():
    for tup in topdicoligs[k]: 
        topligslt_.append('_'.join([tup[0],k])) # ex: 'Ptn_FAPs'
        
print(f"'\n'Getting the receptors for selected topligs list of tuples (topligslt_)") 
deco = {} 
for day in days:
    daydf = hom[day].frame 
    topligdfrec = daydf[daydf['uniq_Ligand_symbol'].isin(topligslt_)]
    #print(topligdfrec[['uniq_Receptor_symbol', 'Edge average expression weight', 'uniq_Ligand_symbol']])
    deco = filldaywvectorstodict(topligdfrec, deco, day)         

print(f"\n Selecting top {NtoprecBylig} receptors for each of the selected ligands")
defdicofiltered = sortandfilterltupdic(reorderdico_fillvariancew(deco), NtoprecBylig, True)


######### find only for sCs the receptors for selected ligands
decosatelli = {}
for day in days:
    daydf = hom[day].frame
    prepadf = daydf[daydf['Target cluster'] == 'sCs']
    topligdfrec = prepadf[prepadf['uniq_Ligand_symbol'].isin(topligslt_)]
    decosatelli = filldaywvectorstodict(topligdfrec, decosatelli,day)
    
defisatellite = sortandfilterltupdic(reorderdico_fillvariancew(decosatelli), 1, True)   


######## do  matrix with signals that arrive to satellite (I start from the end XD)
satellireceptors = set()
ligands = set()

for k in defisatellite.keys():
    lig = k
    ligands.add(k)
    listofreceptors = defisatellite[k]
    for tup in listofreceptors:
        rec = tup[0]
        satellireceptors.add(rec)

# order ligands by celltype origin, alphabetical
preligorder = Sort([i.split("_") for i in ligands], 1, False)
ligorder = ["_".join(l_) for l_ in preligorder]

def createdicoindexes(alistofstrings):
    """
    map of indexes, matches with matrix indexes
    """
    d = {}
    inde = 0
    while inde < len(alistofstrings):
        d[alistofstrings[inde]] = inde
        inde += 1
    return d

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

#ax = sns.heatmap(matrice, cmap="Greens")
#plt.show()
print("we saw that matrix is very very sparse, not the good choice")

"""
better alternative:
    heatmap with weights across time
"""
#print(ligorder)
#print(satellireceptors)
rowlabels = [] 
origcelltype = []
targetcelltype = [] 
day0 = []
day2 = []
day4 = []
day7 = []

for k in ligorder:
    tmpli = k.split("_")
    lig, origct = tmpli[0], tmpli[1]
    print(lig)
    listofreceptors = defisatellite[k] # use key directly to catch lig infos
    for tup in listofreceptors:  # tup
        tmpre = tup[0]
        tmpre = tmpre.split("_")
        rec, targrec = tmpre[0], tmpre[1]
        rowlabels.append(" --> ".join([lig,rec]))
        origcelltype.append(origct)
        targetcelltype.append(targrec)
        day0.append(tup[2][0])
        day2.append(tup[2][1])
        day4.append(tup[2][2])
        day7.append(tup[2][3])
        
#'Index': [i for i in range(len(rowlabels))],
dfplot = pd.DataFrame({  'LR_pairs': rowlabels,
                       '0' : np.log10(day0),
                        '2' : np.log10(day2),
                       '4' : np.log10(day4),
                       '7' : np.log10(day7)})
dfplot = dfplot.set_index(['LR_pairs'])

# celltypes labels:
#row_co = [ (.2, .2,.2, .5) for i in range(84)]
colordictio = {'ECs': '#0000FF', 'M1':'#0000FF',
             'M2':'#0000FF', 'FAPs' : '#0000FF', 'Neutro': '#0000FF', 'sCs': '#0000FF'}
row_co = pd.DataFrame({"Sender_Type" : ['#0000FF' for i in range(84)],
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
for label in colordictio.keys():
    print(label)
    g.ax_col_dendrogram.bar(1.5, 0 , color=colordictio[label],
                            label=label, linewidth=0)
g.ax_col_dendrogram.legend(loc="center right", ncol=2, title="Sender cell types")
ax = g.ax_heatmap
ax.set_xlabel("DAYS")
ax.tick_params( bottom=False) 
#g.cax.set_position([.03, .2, .03, .45])
g.fig.suptitle("Ligand-Receptor pairs: \n top signals to Satellite Cell receptors" ,
                   fontsize=24, horizontalalignment='center' )

g.savefig("MainSignalstoSatelliteCells.pdf", figsize=(12,18))


#g.ax_col_dendrogram.legend(loc="center", ncol=6)
#plt.show()

#g.cax.set_position([.15, .2, .03, .45])
#plt.setp(g.ax_heatmap.xaxis.get_majorticklabels(), fontsize=8)


# import seaborn as sns; sns.set(color_codes=True)
# import string
# iris = sns.load_dataset("iris")
# species = iris.pop("species")
# lut = dict(zip(species.unique(), "rbg"))
# samples = np.repeat(list(string.ascii_letters[0:8]),20)[:150]
# sample_cols = dict(zip(set(samples), sns.color_palette("cubehelix", 8)))

# row_colors = pd.DataFrame({'sample':['g' for i in samples]})
# g = sns.clustermap(iris, row_colors=row_colors,row_cluster=False)

# for k in defisatellite.keys():satellireceptors = []ligands = []


#     lig = k
#     listofreceptors = defisatellite[k]
#     for tup in listofreceptors:
#         rec = tup[0]
#         weis_ = tup[2]
#         delta1 = weis_[1] - weis_[0]
#         delta2 = weis_[2] - weis_[1]
#         delta3 = weis_[3] - weis_[2]
#         satellireceptors.append(rec)
#         ligands.append(lig)
#         delta2vs0.append(delta1)
#         delta4vs2.append(delta2)
#         delta7vs4.append(delta3)
 
# **       
# iris = sns.load_dataset("iris").sample(n=20)
# species = iris.pop("species")
# lut = dict(zip(species.unique(), "rbg"))
# row_colors = species.map(lut)
# g = sns.clustermap(iris, row_colors=row_colors,row_cluster=False)

# df = pd.DataFrame({'Idx1': ['Bar', 'Bar', 'Foo', 'Foo', 'Hop', 'Hop', 'Hop'],
#                     'Idx2': ['a', 'b', 'c', 'd', 'e', 'f', 'g'],
#                     'Col1': np.random.rand(7),
#                     'Col2': np.random.rand(7)})
# df = df.set_index(['Idx1', 'Idx2'])

# g = sns.clustermap(df, center=1, row_cluster=False, cmap="GnBu", yticklabels=True, xticklabels=True, linewidths=0.004)
# g.ax_heatmap.yaxis.set_ticks_position("left")

# plt.setp(g.ax_heatmap.xaxis.get_majorticklabels(), fontsize=12)
# plt.setp(g.ax_heatmap.yaxis.get_majorticklabels(), fontsize=12)
# plt.show()

# **


#cutvar, cutvalue = seevar_setcutWeights(deco, 0.75,0.95)            
  
# def seevar_setcutWeights(dex, QUANT4var, QUANT4valu ):   
#     """ I know i am repeating myself: this is similar to seevariances_setcutoff
#     helps to seek for best cutoffs  ): 
#         - a variance cutoff when more than 2 real products values (not -1)
#         - an product cutoff if only 1 real product available
#     """
#     varthroughdays = [] # variances (if 2 or more day measured)
#     pAlone = [] # single day metrics (no metrics at 2 or more days)
#     for tupk in dex.keys():
#         realvalues = [i for i in dex[tupk] if i > -1]
#         if len(realvalues) > 1:
#             print(realvalues)
#             varthroughdays.append(np.var(realvalues))  
#         elif len(realvalues) == 1 :
#             pAlone.append(realvalues[0])
#     nonzerovars = [i for i in varthroughdays if i > 0.0]
#     autoCutoffvar = np.quantile(np.log10(nonzerovars), q=QUANT4var)    
#     plt.hist(np.log10(nonzerovars), alpha=0.8, color="orangered")
#     plt.axvline(x=autoCutoffvar, color='gray', linestyle='--')
#     plt.title("variances of ... across days")     
#     plt.xlabel("log10(variance)")
#     plt.ylabel("n")
#     plt.show()
#     nonzerop = [ i for i in pAlone if i > 0 ]
#     autoCutoffpAlone = np.quantile(np.log10(nonzerop), q=QUANT4valu)
#     plt.hist(np.log(pAlone), color='gold', alpha=0.6)
#     plt.axvline(x=autoCutoffpAlone, color='orangered', linestyle='--')
#     plt.title(";;;;;; when 1 single day available")
#     plt.xlabel("log10(value)"); plt.ylabel("n")
#     plt.show()
#     return 10**autoCutoffvar, 10**autoCutoffpAlone
  


