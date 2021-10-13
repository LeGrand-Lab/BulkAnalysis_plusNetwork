#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
see header step1 

@author: johaGL
"""
import os
import pickle
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from scipy import stats
import copy
from boxFunClassVersion3 import *

     
def Sort(tups_, myreverse=None):
    if myreverse is None:
        myreverse = False
    tups_.sort(key = lambda x : x[1], reverse=myreverse)
    return(tups_)

def sortandfilterltupdic(groupedbct, N):
    finalgrouped = {}
    for k in groupedbct.keys():    
        givesorted = Sort(groupedbct[k], myreverse=True)
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

def yieldtopUnique(fileligs, allcelltypes, N=None):
    if N is None:
        N = 25
    with open(fileligs, 'r') as f:
        text_ = f.readlines()
    uniquerelsdico = makeUniqueRel(text_)
    dicobyct = groupsdicoUnique(uniquerelsdico, allcelltypes)
    top = sortandfilterltupdic(dicobyct, N)
    for k in top.keys():    
        print(f'  *   {k} top ligands extracted: {len(top[k])}  *  ')
    print()
    return top


lr = pickle.load( open( '../graphobjs/dictio_lr.p', 'rb') )
hom = lr['Young']
print(hom['D0'].frame.columns)
allcelltypes = ['ECs', 'FAPs','M1','M2', 'Neutro', 'sCs']

#fileligs = "selectedLigands.txt"
fileligs = "selectedLigandsVersion3.txt"

topdicoligs = yieldtopUnique(fileligs, allcelltypes, 10)

# get the receptors for selected dico:
# put all L-R  in tuples list: ... put weight  
weights_ = [-1,-1,-1,-1]
daydf = hom['D0'].frame    
ex = daydf.groupby(['Ligand symbol','Sending cluster'])
"""
outputs : dico { (ligand,recep) : [w1,w2,w3,w4]}

"""




