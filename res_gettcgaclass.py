#!/usr/bin/env python
# -*- coding: utf-8 -*-
# @Author : Ji Hongchen
# @Email : jhca.xyt@163.com
# @Last version : 2022-01-03 16:33
# @Filename : res_get_tcga_class.py

import pandas as pd
import numpy as np
import cfg

item = 'TYPE2'
COEF_N = 1
CFG = cfg.cfg()
CLIC_INFO = pd.read_csv('~/Documents/MANU/BI/PANC_ECM/Manuscript/Data/TCGA_ClinicalwithCluster.csv',
                        index_col=[0])

# thresholdlist = [0.483344879565618, 0.484810038619059]  #type1
thresholdlist = [0.52765,
                 0.434, -38.75932283374597] #type

EXP_FILE = pd.read_csv(
    '~/Documents/MANU/BI/PANC_ECM/DATA/SELECTED_MRNA_SYMBOL.csv', index_col=[0])
EXP_FILE = pd.DataFrame(np.array(EXP_FILE).T,
                        index=EXP_FILE.columns, columns=EXP_FILE.index)
EXP_FILE = np.log2(EXP_FILE)

COEF_V = thresholdlist[COEF_N - 1]

COEF = pd.read_csv('~/Documents/MANU/BI/PANC_ECM/Manuscript/Data/' + item + '_lasso_coef' + str(COEF_N) + '.csv',
                   index_col=[0])
print(COEF)
INCP = COEF['s1'].loc['(Intercept)']
COEF_LIST = list(COEF.index)
COEF_LIST.remove('(Intercept)')
COEF_COEF = np.array(COEF.loc[COEF_LIST]).reshape(-1, 1)

EXP_FILE = EXP_FILE[COEF_LIST]

lassoresult = np.array(EXP_FILE) @ COEF_COEF + INCP
EXP_FILE['lasso_pred'] = lassoresult
EXP_FILE['lasso_result'] = 'a'
EXP_FILE['lasso_result'][EXP_FILE['lasso_pred'] >= COEF_V] = 'STROMA-HIGH'
EXP_FILE['lasso_result'][EXP_FILE['lasso_pred'] < COEF_V] = 'STROMA-LOW'
EXP_FILE = EXP_FILE[['lasso_pred', 'lasso_result']]

clinicalinfo = CLIC_INFO.join(EXP_FILE)
print(clinicalinfo)
print('high: ', len(
    clinicalinfo[clinicalinfo['lasso_result'] == 'STROMA-HIGH']))
print('low: ', len(clinicalinfo[clinicalinfo['lasso_result'] == 'STROMA-LOW']))

clinicalinfo.to_csv(
    '~/Documents/MANU/BI/PANC_ECM/Manuscript/Data/Clinica_TCGA_' + item + '_COEF_' + str(COEF_N) + '.csv')
