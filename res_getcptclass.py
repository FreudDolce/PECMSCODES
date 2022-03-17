#!/usr/bin/env python
# -*- coding: utf-8 -*-
# @Author : Ji Hongchen
# @Email : jhca.xyt@163.com
# @Last version : 2021-12-31 10:13
# @Filename : pl_getclass.py

import pandas as pd
import numpy as np
import cfg

CFG = cfg.cfg()
CPT_INFO = '~/Documents/MANU/BI/PANC_ECM/DATA/CPTAC_MRNA_SYMBOL.csv'
CFG.datapath + 'CPT_clinical.csv'

# thresholdlist = [0.483344879565618, 0.484810038619059]  # type1
thresholdlist = [0.52765,
                 0.434, -38.75932283374597]  # type2
item = 'TYPE2'
COEF_N = 2
COEF_V = thresholdlist[COEF_N - 1]


COEF = pd.read_csv(
    '~/Documents/MANU/BI/PANC_ECM/Manuscript/Data/' +
    item + '_lasso_coef' + str(COEF_N) + '.csv',
    index_col=[0])
COEF.drop(COEF[COEF['s1'] == 0].index, inplace=True)
INCP = COEF['s1'].loc['(Intercept)']
COEF_LIST = list(COEF.index)
COEF_LIST.remove('(Intercept)')
COEF_COEF = np.array(COEF.loc[COEF_LIST]).reshape(-1, 1)

REF_EXP = pd.read_csv(
    '~/Documents/MANU/BI/PANC_ECM/DATA/SELECTED_MRNA_SYMBOL.csv',
    index_col=[0])
REF_EXP = pd.DataFrame(np.array(REF_EXP).T,
                       index=REF_EXP.columns, columns=REF_EXP.index)
REF = REF_EXP[COEF_LIST]
REF = np.log2(REF)
col_ref_m = np.array(REF).mean()
print(col_ref_m)
col_ref_s = np.array(REF).std()

CPT_CLINIC = pd.read_csv(
    '~/Documents/MANU/BI/PANC_ECM/Manuscript/Data/CPT_clinical.csv',
    index_col=[0])

if __name__ == '__main__':
    cpt_info = pd.read_csv(CFG.CPT_EXP_FILE, index_col=[0])
    # cpt_info = pd.DataFrame(np.array(cpt_info).T,
    #                       index=cpt_info.columns,
    #                       columns=cpt_info.index)
    cpt_info = cpt_info[COEF_LIST]
    cpt_info = np.log2(cpt_info)
    cpt_col_m = np.array(cpt_info).mean()
    cpt_col_s = np.array(cpt_info).std()
    cpt_info = (cpt_info - cpt_col_m) / cpt_col_s
    #cpt_info = cpt_info * cpt_col_s + cpt_col_m
    cpt_info = cpt_info * col_ref_s + col_ref_m
    print(np.array(cpt_info).mean())
    cpt_result = np.array(
        cpt_info) @ np.array(COEF_COEF).reshape(-1, 1) + INCP
    cpt_info['lasso_pred'] = cpt_result
    cpt_info['lasso_result'] = 'a'
    cpt_info['lasso_result'][cpt_info['lasso_pred'] >= COEF_V] = 'STROMA-HIGH'
    cpt_info['lasso_result'][cpt_info['lasso_pred'] < COEF_V] = 'STROMA-LOW'
    cpt_info = cpt_info[['lasso_pred', 'lasso_result']]
    print('high: ', len(cpt_info[cpt_info['lasso_result'] == 'STROMA-HIGH']))
    print('low: ', len(cpt_info[cpt_info['lasso_result'] == 'STROMA-LOW']))

    cpt_with_clic = CPT_CLINIC.join(cpt_info)
    cpt_with_clic.to_csv(
        '~/Documents/MANU/BI/PANC_ECM/Manuscript/Data/Clinica_CPT_' + item + '_COEF_' + str(COEF_N) + '.csv')
