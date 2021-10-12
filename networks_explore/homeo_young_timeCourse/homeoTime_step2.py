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
from boxFunClass import *


lr = pickle.load( open( '../graphobjs/dictio_lr.p', 'rb') )

with open("selectedLigands.txt", 'r') as f:
    text_ = f.readlines()


allcelltypes = ['ECs', 'FAPs','M1','M2', 'Neutro', 'sCs']
groupedbct = {}
finalgrouped = {}
for i in allcelltypes:
    groupedbct[i] = []
    finalgrouped[i] = []
    
print(text_[0]) # the header: symbol	cellty	varianceORexpression	criterionselection
for i in range(1,len(text_)):
    tutu = tuple(text_[i].strip().split('\t'))
    groupedbct[tutu[1]].append((tutu[0], float(tutu[2])))

def Sort(tups_, myreverse=None):
    if myreverse is None:
        myreverse = False
    tups_.sort(key = lambda x : x[1], reverse=myreverse)
    return(tups_)

# parse this dico of tuples, rank by second element, choose first N :  
N = 25
for k in groupedbct.keys():    
    print(f'  *   {k} top ligands are: {len(groupedbct[k])}  *  ')
    print(type(groupedbct[k][-1][1]))
    givesorted = Sort(groupedbct[k], myreverse=True)
    if len(givesorted) > N:
        finalgrouped[k] = givesorted[:N]
    else:
        finalgrouped[k] = givesorted

    
"""
separate into lists by celltype:
  *   Neutro top ligands are: 28  *  
[('Icam1', 266.706), ('Il1b', 4022.527), ('Cd14', 21...   ]
None
  *   sCs top ligands are: 22  *  
[('Serping1', 1157.37), ('Col1a2', 3402.00), ('Igf2', 585 ... ]                                      
"""

 



