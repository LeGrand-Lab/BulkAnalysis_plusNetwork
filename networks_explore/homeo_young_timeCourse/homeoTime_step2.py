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
from boxFunClassVersion3 import *

############## find for all celltypes the receptors for selected ligands

lr = pickle.load( open( '../graphobjs/dictio_lr.p', 'rb') )
hom = lr['Young']
days = ['D0', 'D2', 'D4', 'D7']
allcelltypes = ['ECs', 'FAPs','M1','M2', 'Neutro', 'sCs']

#fileligs = "selectedLigands.txt"
fileligs = "selectedLigandsVersion3.txt"

topdicoligs = yieldtopUnique(fileligs, allcelltypes, 10)
toplol_ = []
for k in topdicoligs.keys():
    for tup in topdicoligs[k]: 
        toplol_.append('_'.join([tup[0],k])) # ex: 'Ptn_FAPs'
        
# get the receptors for selected dico:
# put all L-R  in tuples list: ... put weight  
def banana(sue, deco):
    weis_ = [-1,-1,-1,-1]    
    dki = {'D0':0,'D2':1, 'D4':2, 'D7':3}  
    for i, row in sue.iterrows():
        fromvertex = row['uniq_Ligand_symbol']
        tovertex = row['uniq_Receptor_symbol']
        weiday = row['Edge average expression weight']
        try:
            deco[(fromvertex,tovertex)][dki[day]] = float(weiday)
        except KeyError:
            deco[(fromvertex,tovertex)] = copy.deepcopy(weis_)
            deco[(fromvertex,tovertex)][dki[day]] = float(weiday)  
    return deco

    
def pineaple(deco):
    varsdico = {}
    for tupk in deco.keys():
        """
        output: dico { ligand : [(rec1, varianceweight),(rec2, varianceweight),... }
        """
        ligand, receptor = tupk[0], tupk[1]
        truewei = [i for i in deco[tupk] if i > -1]
        if len(truewei) > 1:
            v = np.var(truewei)
        elif len(truewei) == 1:
            v = truewei[0]
        try :
            varsdico[ligand].append((receptor, v ))
        except KeyError:
            varsdico[ligand] = [(receptor, v)]
    return varsdico
 
deco = {} 
for day in days:
    daydf = hom[day].frame 
    """
    output : dico { (ligand,recep) : [w1,w2,w3,w4]}
    
    """
    sue = daydf[daydf['uniq_Ligand_symbol'].isin(toplol_)]
    #print(sue[['uniq_Receptor_symbol', 'Edge average expression weight', 'uniq_Ligand_symbol']])
    deco = banana(sue, deco)         

defdicofiltered = sortandfilterltupdic(pineaple(deco), 2, True)

######### find only for sCs the receptors for selected ligands
decosatelli = {}
for day in days:
    daydf = hom[day].frame
    presue = daydf[daydf['Target cluster'] == 'sCs']
    sue = presue[presue['uniq_Ligand_symbol'].isin(toplol_)]
    decosatelli = banana(sue, decosatelli)
    
defisatellite = sortandfilterltupdic(pineaple(decosatelli), 3, True)   


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
  


