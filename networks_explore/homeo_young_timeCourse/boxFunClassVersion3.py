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

# used in step 1 and 2

class LRinfos:  
    """
    class to handle  Natmi dataframe results
    one object by result ! 
    to create object, predf is needed (the Edges opened with pandas csv)
    use 'frame' attribute to get dataframe (is a pandas object)
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
    EXPLANATION: 4 float list (1 expresion*specificity value per day)
            An artifactual -1 means no data for that symbol+celltype+day
            """
    dki = {'D0':0,'D2':1, 'D4':2, 'D7':3}  
    foo = [-1,-1,-1,-1]     
    ex = daydf.groupby(['Ligand symbol','Sending cluster'])
    ex = ex[['Ligand symbol','Sending cluster',   
             'Ligand average expression value',           
             'Ligand derived specificity of average expression value' ]]
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
         PRODUCTO = row['Ligand derived specificity of average expression value'] *\
             row['Ligand average expression value'] 
         dex[sym][row['Sending cluster']][dki[day]] = PRODUCTO
    return dex
                
def seevariances_setcutoff(dex, QUANT4var, QUANT4prod ):   
    """
    helps to seek for best cutoffs : 
        - a variance cutoff when more than 2 real products values (not -1)
        - an product cutoff if only 1 real product available
    """
    varthroughdays = [] # variances (if 2 or more day measured)
    pAlone = [] # single day metrics (no metrics at 2 or more days)
    for k in dex.keys():
        for ct in dex[k]:
            realvalues = [float(i) for i in dex[k][ct] if i > -1]
            if len(realvalues) > 1:
                print(realvalues)
                varthroughdays.append(np.var(realvalues))  
            elif len(realvalues) == 1 :
                pAlone.append(realvalues[0])
    nonzerovars = [i for i in varthroughdays if i > 0.0]
    autoCutoffvar = np.quantile(np.log10(nonzerovars), q=QUANT4var)    
    plt.hist(np.log10(nonzerovars), alpha=0.6)
    plt.axvline(x=autoCutoffvar, color='gray', linestyle='--')
    plt.title("variances of Ligands' (expression*specificities) across days")     
    plt.xlabel("log10(variance)")
    plt.ylabel("n")
    plt.show()
    nonzerop = [ i for i in pAlone if i > 0 ]
    autoCutoffpAlone = np.quantile(np.log10(nonzerop), q=QUANT4prod)
    plt.hist(np.log(pAlone), color='green', alpha=0.4)
    plt.axvline(x=autoCutoffpAlone, color='gray', linestyle='--')
    plt.title("Expression*specificities when 1 single day available")
    plt.xlabel("log10(value)"); plt.ylabel("n")
    plt.show()
    return 10**autoCutoffvar, 10**autoCutoffpAlone

def pickovercutoff(dex, varcutoff, expcutoff):
    """
    yields a list of selected (over cutoffs) tuples : 
        (a_symbol, a_celltype, the_varianceORmetric, a_criterion_cutoff)
    """
    selection = []
    for k in dex.keys():
        for ct in dex[k]:
            realvalues = [i for i in dex[k][ct] if i > -1]
            if len(realvalues) > 1 :
                tmpvariance = np.var(realvalues)
                if tmpvariance > varcutoff:
                    selection.append((k, ct, tmpvariance, "VARIANCEBASED"))
            elif len(realvalues) == 1 :
                if realvalues[0] > expcutoff:
                    selection.append((k, ct, realvalues[0], "prodbased"))
    return selection

# used in step 2 :
def Sort(tups_, pos,  myreverse=None):
    if myreverse is None:
        myreverse = False
    tups_.sort(key = lambda x : x[pos], reverse=myreverse)
    return(tups_)

def sortandfilterltupdic(groupedbct, N, myreverse=None):
    """"
    input : {'FAPs' : [('gene1', value), ('gene2',value) ... ] , .... }
    where values are variances
    by default sorts from bigger to smaller tuple by 2nd element values (True)
    Output: same dictionnary but N (number) selected (gene,value) tuples
    """
    if myreverse is None:
        myreverse = True
    finalgrouped = {}
    for k in groupedbct.keys():    
        givesorted = Sort(groupedbct[k], 1, myreverse) 
        if len(givesorted) > N:
            finalgrouped[k] = givesorted[:N]
        else: # fewer than N tuples
            finalgrouped[k] = givesorted
    return finalgrouped

def makeUniqueRel(text_):
    """
    when seing .txt : multiple cases of same symbol repeated on several celltypes: 
        Psap	ECs	41461.86834053416	VARIANCEBASED
    Psap	FAPs	347.32303895593446	VARIANCEBASED
    Psap	M1	44467.65913367455	VARIANCEBASED
    Psap	M2	27242.193053371284	VARIANCEBASED    
    solution: pick max celltype value for a given gene symbol (uniqueness)
    """
    symbolstocellty = {}
    print(text_[0]) # the header: symbol	cellty	...	criterionselection
    for i in range(1,len(text_)):
        tutu = tuple(text_[i].strip().split('\t'))
        try:
            symbolstocellty[tutu[0]].append((tutu[1],tutu[2]))
        except KeyError:
            symbolstocellty[tutu[0]] = [(tutu[1],tutu[2])]
    themaxdico = sortandfilterltupdic(symbolstocellty, 1)
    return themaxdico

def groupsdicoUnique(uniquerelsdico, allcelltypes):
    """
    organize dico by gene, into dico by celltype:
       {  'Neutro' : 
        [('Icam1', 266.706), ('Il1b', 4022.527), ('Cd14', 21...   ] , 
         'sCs'  :
        [('Serping1', 1157.37), ('Col1a2', 3402.00), ('Igf2', 585 ... ]  ... }                                    
    """
    groupedbct = {}
    for i in allcelltypes:
        groupedbct[i] = []
    for k in uniquerelsdico.keys():
        tup = uniquerelsdico[k][0] # [(tuple)] take only tuple
        groupedbct[tup[0]].append((k, tup[1]))
    for k in groupedbct.keys():    
        print(f'  *   {k} top ligands found: {len(groupedbct[k])}  *  ')
    print()
    return groupedbct

def yieldtopUnique(fileligs, allcelltypes, N=None):
    """"outputs a dictionary with keys celltypes:
        { 'sCs' : [('Vim', '9538.6709'), ('Jam3', '93.680415'),  ...   
    """
    if N is None:
        N = 25
    with open(fileligs, 'r') as f:
        text_ = f.readlines()
    uniquerelsdico = makeUniqueRel(text_)
    dicobyct = groupsdicoUnique(uniquerelsdico, allcelltypes)
    top = sortandfilterltupdic(dicobyct, N)
    for k in top.keys():    
        print(f'  *   {k} top ligands extracted: {len(top[k])}  *  ')    
    return top

def seeNonUnique_andgroup(text_, allcelltypes):
    groupedbct = {}
    for i in allcelltypes:
        groupedbct[i] = []
    print(text_[0]) # the header: symbol	cellty	...	criterionselection
    for i in range(1,len(text_)):
        tutu = tuple(text_[i].strip().split('\t'))
        groupedbct[tutu[1]].append((tutu[0], float(tutu[2])))  
    for k in groupedbct.keys():    
        print(f'  *   {k} top ligands are: {len(groupedbct[k])}  *  ')
    print()
    return groupedbct

def yieldtopNonUnique(fileligs, allcelltypes, N=None):
    if N is None:
        N = 25
    with open(fileligs, 'r') as f:
        text_ = f.readlines()
    groupedltup = seeNonUnique_andgroup(text_, allcelltypes)
    final = sortandfilterltupdic(groupedltup, N)
    return final








