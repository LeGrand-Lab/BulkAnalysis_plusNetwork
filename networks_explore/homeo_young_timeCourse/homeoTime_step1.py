#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""

I want to extract L-R pairs which, accross time, change the most.
 BUT celltypes are not the same accross days !
 how to handle ?  : 
     step 1 : 
         - create sender NESTED dict using Natmi output dataframes all days :  
           external keys are the gene symbols
            calculate the product expression*specificity, nested dictionary,    
                  keys are celltypes, values are list product values (4values)
                  calculate variance of generated values across days
        - select, across sequential days,  MOST VARIABLE ligands (variance cutoff)
              (Note: if value for only one day, filter by simple value cutoff)
     step 2 :
         - see homeoTime_step2.py
     
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

hom = lr['Young']
days = ['D0', 'D2', 'D4', 'D7']
allcelltypes = [ "ECs" ,
    "FAPs"  ,
    "Inflammatory-Mac"  ,
    "Resolving-Mac"  ,
    "Neutrophils" ,
    "MuSCs" ]
ds = {}
for k in days:
    df1 = hom[k].frame  # hom[k].frame.sample(n=100) helps tests
    ds = filldicoSenders(ds, df1, allcelltypes,  k)  # yes, add days iteratively to same dictionnary
    
varcutoff, expcutoff = seevariances_setcutoff(ds, 0.65, 0.75) # QUANTILEs best tuned
print(varcutoff)
mysele = pickovercutoff(ds, varcutoff, expcutoff)
texto = "symbol\tcellty\tvarianceORproduct\tcriterionselection\n"
for i in mysele:
    texto += '\t'.join([str(x) for x in i]) 
    texto += '\n'               

with open("selectedLigandsVersion3.txt", "w") as f:
    f.write(texto)
    f.close()
    











 



