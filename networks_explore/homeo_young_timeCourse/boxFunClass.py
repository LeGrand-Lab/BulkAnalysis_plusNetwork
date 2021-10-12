#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
functions necessary to step1, 2, etc
for Homeostasis Time course L-R analysis

@author: johaGL, oct 2021
"""

import os
import pickle
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from scipy import stats
import copy


class LRinfos:  
    """
    class to handle  Natmi dataframe results
    one object by result ! 
    to create object, predf is needed (the Edges opened with pandas csv)
    use 'frame' attribute to get dataframe suitable for graph conversion
    """
    def __init__(self, age, day, predf):
        self.age = age
        self.day = day
        self.predf  = predf
        self.makeunique_symbo_cellty()
        self.about = f"object age {age}, day {day}, use 'frame' attribute for more!"
        
    def makeunique_symbo_cellty(self):
        otab = self.predf
        otab['uniq_Ligand_symbol'] = otab['Ligand symbol'] + '_' + otab['Sending cluster']
        otab['uniq_Receptor_symbol'] = otab['Receptor symbol'] + '_' + otab['Target cluster']
        self.frame = otab  # this adds attribute 'frame'     
        
    def filterZero(self):  
        tmp = self.frame.loc[self.frame['Edge average expression derived specificity'] > 0]
        self.frame = tmp  # yield only non zero edges dataframe
    
    def filterOnEdgeslog10(self, cutoff):
        if min(self.frame['Edge average expression derived specificity']) < 0:
            self.filterZero()
        self.frame['log10_edge_sp'] = np.log10(np.array(self.frame['Edge average expression derived specificity']))
        tf = self.frame.loc[ self.frame['log10_edge_sp'] >= cutoff ]
        self.filtered = tf


def filldicoSenders(dex, daydf, allcelltypes, day): 
    """
    Parameters
    ----------
    dex : dictionary
    daydf : pandas L-R dataframe
    allcelltypes: all possible cell types
    day : string (ex. 'D7')
        
    Returns
    -------
    dex : dicionary 
       'Sostdc1': {'ECs': [4.12, 0.39, 2.30, 40.75], 'FAPs': [0. ...],
    EXPLANATION: 4 float list (1 expresion value per day)
            An artifactual -1 means no data for that symbol+celltype+day
            """
    dki = {'D0':0,'D2':1, 'D4':2, 'D7':3}  
    foo = [-1,-1,-1,-1]     
    ex = daydf.groupby(['Ligand symbol','Sending cluster'])
    ex = ex[['Ligand symbol','Sending cluster','Ligand average expression value']]
    ex = ex.apply(pd.DataFrame)
    ex.index = range(ex.shape[0])   
    allsymbols = list(pd.unique(ex['Ligand symbol']))
    nestedd = {}
    for t in allcelltypes:
        nestedd[t] = copy.deepcopy(foo)
    for s in allsymbols:
        if s not in dex.keys():
            dex[s] = copy.deepcopy(nestedd)
    for i,row in ex.iterrows():
         print(row)
         sym = str(row['Ligand symbol']) 
         dex[sym][row['Sending cluster']][dki[day]] = row['Ligand average expression value']
    return dex
                
def seevariances_setcutoff(dex, QUANT4var, QUANT4exp ):   
    """
    helps to seek for best cutoffs : 
        - a variance cutoff when more than 2 real expression values (not -1)
        - an expression cutoff if only 1 real expression
    """
    varthroughdays = [] # variances (if 2 or more day measured)
    expressAlone = [] # single day expression (no metrics at 2 or more days)
    for k in dex.keys():
        for ct in dex[k]:
            realvalues = [i for i in dex[k][ct] if i > -1]
            if len(realvalues) > 1:
                print(realvalues)
                varthroughdays.append(np.var(realvalues))  
            elif len(realvalues) == 1 :
                expressAlone.append(realvalues)
    nonzerovars = [i for i in varthroughdays if i > 0]
    autoCutoffvar = np.quantile(np.log10(nonzerovars), q=QUANT4var)    
    plt.hist(np.log10(nonzerovars), alpha=0.6)
    plt.axvline(x=autoCutoffvar, color='gray', linestyle='--')
    plt.title("variances of Ligands expression across days")     
    plt.xlabel("variance")
    plt.ylabel("n")
    plt.show()
    autoCutoffExpAlone = np.quantile(np.log10(expressAlone), q=QUANT4exp)
    plt.hist(np.log(expressAlone), color='green', alpha=0.4)
    plt.axvline(x=autoCutoffExpAlone, color='gray', linestyle='--')
    plt.show()
    return 10**autoCutoffvar, 10**autoCutoffExpAlone

def pickovercutoff(dex, varcutoff, expcutoff):
    """
    yields a list of selected (over cutoffs) tuples : 
        (a_symbol, a_celltype, the_varianceORExpression, a_criterion_cutoff)
    """
    selection = []
    for k in dex.keys():
        for ct in dex[k]:
            realvalues = [i for i in dex[k][ct] if i > -1]
            if len(realvalues) > 1 :
                tmpvariance = np.var(realvalues)
                if tmpvariance > varcutoff:
                    selection.append((k, ct, tmpvariance, "variancebased"))
            elif len(realvalues) == 1 :
                if realvalues[0] > expcutoff:
                    selection.append((k, ct, realvalues[0], "EXPBASED"))
    return selection










