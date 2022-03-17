#!/usr/bin/env python
# -*- coding: utf-8 -*-
# @Author : Ji Hongchen
# @Email : jhca.xyt@163.com
# @Last version : 2022-03-07 18:00
# @Filename : res_devide_scclinical_data.py

import pandas as pd
import numpy as np
import os

GENELIST = ['COL17A1', 'AREG', 'KLHL32',
            'CDA', 'POSTN', 'SLC2A1', 'FN1', 'INHBA']
REF_RNA = pd.read_csv('~/Documents/MANU/BI/PANC_ECM/DATA/SELECTED_MRNA_SYMBOL.csv',
                      index_col=[0])
REF_RNA = pd.DataFrame(np.array(REF_RNA).T,
                       index=REF_RNA.columns,
                       columns=REF_RNA.index)
REF_RNA = REF_RNA[GENELIST]
REF_RNA = np.log2(REF_RNA)

clinical_exp = pd.read_csv('~/Documents/MANU/BI/PANC_ECM/DATA/Clinical_data_exp/gene_expression.csv',
                           index_col=[0])

clinical_exp = clinical_exp[GENELIST]
clinical_exp = np.log2(clinical_exp)
clinical_exp.replace(-np.inf, np.nan, inplace=True)
clinical_exp.fillna(0, inplace=True)

for g in GENELIST:
    ref_mean = np.array(REF_RNA[g]).mean()
    ref_std = np.array(REF_RNA[g]).std()
    clic_mean = np.array(clinical_exp[g]).mean()
    clic_std = np.array(clinical_exp[g]).std()

    clinical_exp[g] = (clinical_exp[g] - clic_mean) / clic_std
    clinical_exp[g] = clinical_exp[g] * ref_std + ref_mean

COEF = pd.read_csv('~/Documents/MANU/BI/PANC_ECM/Manuscript/Data/TYPE2_lasso_coef2.csv',
                   index_col=[0])
COEF.drop(COEF[COEF['s1'] == 0].index, inplace=True)
INCP = COEF['s1'].loc['(Intercept)']
COEF_LIST = list(COEF.index)
COEF_LIST.remove('(Intercept)')
COEF_COEF = np.array(COEF.loc[COEF_LIST]).reshape(-1, 1)

clic_pred = np.array(clinical_exp) @ np.array(COEF_COEF).reshape(-1, 1) + INCP
clinical_exp['pred'] = clic_pred
clinical_exp['class'] = 0
clinical_exp['class'][clinical_exp['pred'] > 0.434] = 1
print(clinical_exp)
clinical_exp.to_csv(
    '~/Documents/MANU/BI/PANC_ECM/DATA/Clinical_data_exp/gene_expression.csv')
