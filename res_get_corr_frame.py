#!/usr/bin/env python
# -*- coding: utf-8 -*-
# @Author : Ji Hongchen
# @Email : jhca.xyt@163.com
# @Last version : 2022-01-20 15:47
# @Filename : res_get_corr_frame.py

import pandas as pd
import numpy as np
import os

EXP = pd.read_csv('~/Documents/MANU/BI/PANC_ECM/DATA/SELECTED_MRNA_SYMBOL.csv',
                  index_col=[0])
EXP = np.log2(EXP)
for c in EXP.columns:
    exp_m_c = EXP[c].mean()
    exp_m_s = EXP[c].std()
    EXP[c] = (EXP[c] - exp_m_c) / exp_m_s
EXP = pd.DataFrame(np.array(EXP).T, index=EXP.columns, columns=EXP.index)
genelist = list(pd.read_csv('~/Documents/MANU/BI/PANC_ECM/Manuscript/Data/TYPE2_lasso_coef2.csv',
                            index_col=[0]).index)
genelist.remove('(Intercept)')

Metalist = list(pd.read_csv('~/Documents/MANU/BI/PANC_ECM/Manuscript/Data/Corr_ros_coef.csv',
                            index_col=[0]).columns)
metadata = EXP[Metalist]

EXP = EXP[genelist]
# Get gsva corr
Gsvaresult = pd.read_csv('~/Documents/MANU/BI/PANC_ECM/DATA/Immune_type_gsva_result.csv',
                         index_col=[0])
gsvalist = list(Gsvaresult.columns)

Gsvaframe = pd.DataFrame(columns=gsvalist, index=genelist)
for c in gsvalist:
    for g in genelist:
        corr_v = np.corrcoef(EXP[g], Gsvaresult[c])
        Gsvaframe[c][g] = corr_v[0][1]

Gsvaframe.to_csv(
    '~/Documents/MANU/BI/PANC_ECM/Manuscript/Data/Draw_corr_gsva.csv')

NU_EXP = pd.read_csv(
    '~/Documents/MANU/BI/PANC_ECM/DATA/SELECTED_MRNA.csv', index_col=[0])
NU_EXP = pd.DataFrame(np.array(NU_EXP).T,
                      index=NU_EXP.columns, columns=NU_EXP.index)
NU_EXP = np.log2(NU_EXP)
NU_EXP.replace(-np.inf, np.nan, inplace=True)
NU_EXP.fillna(0, inplace=True)
for c in NU_EXP:
    nu_c_m = NU_EXP[c].mean()
    nu_c_s = NU_EXP[c].std()
    NU_EXP[c] = (NU_EXP[c] - nu_c_m) / nu_c_s

# Get immune label corr
Drawlist = [5133, 29126, 80380, 1493, 941,
            942, 84868, 201633, 3604, 939, 9760, 953]
im_exp = NU_EXP[Drawlist]
print(im_exp)

Immunelabelframe = pd.DataFrame(columns=Drawlist, index=genelist)
for c in Drawlist:
    for g in genelist:
        immune_c = im_exp[c]
        gene_g = EXP[g]
        corr_v = np.corrcoef(immune_c, gene_g)
        Immunelabelframe[c][g] = corr_v[0][1]

Immunelabelframe.to_csv(
    '~/Documents/MANU/BI/PANC_ECM/Manuscript/Data/Draw_corr_immune.csv')

# Get immune invastion corr
Cibersortresult = pd.read_csv('~/Documents/MANU/BI/PANC_ECM/DATA/CIBERSORT_RESULT_for_plot.csv',
                              index_col=[0])
Drawlist = list(Cibersortresult.columns)
cibersortframe = pd.DataFrame(columns=Drawlist, index=genelist)
for c in Drawlist:
    for g in genelist:
        cibersortframe[c][g] = np.corrcoef(Cibersortresult[c], EXP[g])[0][1]
cibersortframe.to_csv(
    '~/Documents/MANU/BI/PANC_ECM/Manuscript/Data/Draw_corr_ciber.csv')

# Get metabolism frame
metaframe = pd.DataFrame(columns=Metalist, index=genelist)
for c in Metalist:
    for g in genelist:
        metaframe[c][g] = np.corrcoef(metadata[c], EXP[g])[0][1]

metaframe.to_csv(
    '~/Documents/MANU/BI/PANC_ECM/Manuscript/Data/Draw_corr_meta.csv')

# Get drug sensitive frame
Drugdata = pd.read_csv('~/Documents/MANU/BI/PANC_ECM/DATA/DrugSensitive_predict.csv',
                       index_col=[0])
druglist = list(Drugdata.columns)

drugframe = pd.DataFrame(columns=druglist, index=genelist)
for c in druglist:
    for g in genelist:
        drugframe[c][g] = np.corrcoef(Drugdata[c], EXP[g])[0][1]

drugframe.to_csv('~/Documents/MANU/BI/PANC_ECM/Manuscript/Data/Draw_corr_drugsen.csv')
