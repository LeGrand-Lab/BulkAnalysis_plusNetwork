#!/usr/bin/env python3
# -*- coding: utf-8 -*-
import os

NatmiOutDirectory = "/home/mopa/Documents/INMG/BulkAnalysis_plusNetwork/NatmiData/natmiOut_TPM"
Days=["D0","D2","D4","D7"]
TypeCell=["ECs","FAPs","MuSCs","Neutrophils","Inflammatory-Mac","Resolving-Mac"]
intDirectory = "Delta_CltPair_exp_0_spe_0_det_0.3_top_0_signal_lrc2p_weight_mean"
import glob
import time
import pandas as pd
import xlwt
from xlwt.Workbook import *
from pandas import ExcelWriter
import xlsxwriter
import numpy as np


writer = pd.ExcelWriter(NatmiOutDirectory+'/Diff_LR_pop_up_down_regulation_total.xlsx', engine='xlsxwriter')
for D in Days:
    files=glob.glob(NatmiOutDirectory+"/Diff_"+D+"/"+intDirectory+'/*.csv')
    for ct in TypeCell:
        print(ct)
        tableSheetconcat=pd.DataFrame({'Sending cluster':[],'Ligand symbol':[],'Receptor symbol':[],'Target cluster':[],'Delta ligand detection rate':[],'Delta ligand average expression value':[],'Delta ligand derived specificity of average expression value':[],'Delta receptor detection rate':[],'Delta receptor average expression value':[],'Delta receptor derived specificity of average expression value':[],'Edge delta average expression weight':[],'Edge delta average expression derived specificity':[]})
        for f in files:
            motif=ct+"_to"
            if motif in f:
                tableSheet=pd.read_csv(f)
                diffstatus=f.split('/')[-1].split('_')[0] 
                receptor=f.split('/')[-1].split('_')[4]
                if receptor == "to":
                    receptor=f.split('/')[-1].split('_')[5]
                namediffstatus=diffstatus+"_"+ct+"_to_"+receptor
                print(f)
                print(namediffstatus)
                tableSheet["DifferenceStatus"]=np.repeat(namediffstatus,tableSheet.shape[0])
                tableSheetconcat=pd.concat([tableSheetconcat,tableSheet])
        if tableSheetconcat.shape[0] > 2:
            tableSheetconcat.to_excel(writer,sheet_name=D+'_'+ct)
writer.close()
    
"""
NatmiOutDirectory = "/home/bioinfo/BulkAnalysis_plusNetwork/NatmiData/natmiOut_CountNormalised"
Days=["D0","D2","D4","D7"]
TypeCell=["ECs","FAPs","MuSCs","Neutrophils","Inflammatory-Mac","Resolving-Mac"]
intDirectory = "Delta_CltPair_exp_0_spe_0_det_0.3_top_0_signal_lrc2p_weight_mean"

writer = pd.ExcelWriter(NatmiOutDirectory+'/Diff_LR_pop_up_down_regulation_total_CN.xlsx', engine='xlsxwriter')
for D in Days:
    files=glob.glob(NatmiOutDirectory+"/Diff_"+D+"/"+intDirectory+'/*.csv')
    for ct in TypeCell:
        print(ct)
        tableSheetconcat=pd.DataFrame({'Sending cluster':[],'Ligand symbol':[],'Receptor symbol':[],'Target cluster':[],'Delta ligand detection rate':[],'Delta ligand average expression value':[],'Delta ligand derived specificity of average expression value':[],'Delta receptor detection rate':[],'Delta receptor average expression value':[],'Delta receptor derived specificity of average expression value':[],'Edge delta average expression weight':[],'Edge delta average expression derived specificity':[]})
        for f in files:
            motif=ct+"_to"
            if motif in f:
                tableSheet=pd.read_csv(f)
                diffstatus=f.split('/')[-1].split('_')[0] 
                receptor=f.split('/')[-1].split('_')[4]
                if receptor == "to":
                    receptor=f.split('/')[-1].split('_')[5]
                namediffstatus=diffstatus+"_"+ct+"_to_"+receptor
                print(f)
                print(namediffstatus)
                tableSheet["DifferenceStatus"]=np.repeat(namediffstatus,tableSheet.shape[0])
                tableSheetconcat=pd.concat([tableSheetconcat,tableSheet])
        if tableSheetconcat.shape[0] > 2:
            tableSheetconcat.to_excel(writer,sheet_name=D+'_'+ct)
writer.save()

"""

