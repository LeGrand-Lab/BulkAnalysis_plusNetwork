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
     - fill a 3D matrix, a cell in the matrix will be a list of 4 values
                 corresponding to 'Edge average expression weight' by each day
                  (for each day dataframe, get sender-receiver as named in matrix
                   matching to celltype as said in dictionary Z)
        - variance matrix : find the days edges variance from 3D matrix
            - explore  : do some plots of that variance 
            - filter pairs having the biggest variance
            - do plot
        - write a simple csv file containing most variable edges as found 
        in variance matrix, containing same columns as Natmi output dataframes

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

lr = pickle.load( open( '../graphobjs/dictio_lr.p', 'rb') )

hom = lr['Young']
days = ['D0', 'D2', 'D4', 'D7']
allcelltypes = ['ECs', 'FAPs','M1','M2', 'Neutro', 'sCs']
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
    











 



