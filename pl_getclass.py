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
CPT_INFO = CFG.datapath + 'CPT_clinical.csv'
COL_TH = CFG.COL_THERS
HA_TH = CFG.HA_THERS
COL_COEF = pd.read_csv(CFG.lassoresultpath +
                       'SCORE/COL_COEF_2.csv', index_col=[0])
COL_INCP = COL_COEF['s1'].loc['(Intercept)']
COL_COEF_LIST = list(COL_COEF.index)
COL_COEF_LIST.remove('(Intercept)')
COL_COEF_COEF = np.array(COL_COEF.loc[COL_COEF_LIST]).reshape(-1, 1)
HA_COEF = pd.read_csv(CFG.lassoresultpath +
                      'SCORE/HA_COEF_2.csv', index_col=[0])
HA_INCP = HA_COEF['s1'].loc['(Intercept)']
HA_COEF_LIST = list(HA_COEF.index)
HA_COEF_LIST.remove('(Intercept)')
HA_COEF_COEF = np.array(HA_COEF.loc[HA_COEF_LIST]).reshape(-1, 1)

REF_EXP = pd.read_csv(CFG.EXPRESSION_FILE, index_col=[0])
REF_EXP = pd.DataFrame(np.array(REF_EXP).T,
                       index=REF_EXP.columns, columns=REF_EXP.index)
COL_REF = REF_EXP[COL_COEF_LIST]
COL_REF = np.log2(COL_REF)
col_ref_m = np.array(COL_REF).mean()
col_ref_s = np.array(COL_REF).std()
HA_REF = REF_EXP[HA_COEF_LIST]
HA_REF = np.log2(HA_REF)
ha_ref_m = np.array(HA_REF).mean()
ha_ref_s = np.array(HA_REF).std()


if __name__ == '__main__':
    cpt_info = pd.read_csv(CFG.CPT_EXP_FILE, index_col=[0])
    cpt_info = pd.DataFrame(np.array(cpt_info).T,
                            index=cpt_info.columns,
                            columns=cpt_info.index)
    cpt_col_info = cpt_info[COL_COEF_LIST]
    cpt_col_info = np.log2(cpt_col_info)
    cpt_col_m = np.array(cpt_col_info).mean()
    cpt_col_s = np.array(cpt_col_info).std()
    cpt_col_info = (cpt_col_info - cpt_col_m) / cpt_col_s
    cpt_col_info = cpt_col_info * col_ref_s + col_ref_m
    cpt_col_result = np.array(
        cpt_col_info) @ np.array(COL_COEF_COEF).reshape(-1, 1) + COL_INCP
    cpt_col_info['col_pred'] = cpt_col_result
    cpt_col_info['COL'] = cpt_col_result
    cpt_col_info['COL'][cpt_col_info['COL'] >= COL_TH] = 1
    cpt_col_info['COL'][cpt_col_info['COL'] < COL_TH] = 0

    cpt_ha_info = cpt_info[HA_COEF_LIST]
    cpt_ha_info = np.log2(cpt_ha_info)
    cpt_ha_m = np.array(cpt_ha_info).mean()
    cpt_ha_s = np.array(cpt_ha_info).std()
    cpt_ha_info = (cpt_ha_info - cpt_ha_m) / cpt_ha_s
    cpt_ha_info = cpt_ha_info * ha_ref_s + ha_ref_m
    cpt_ha_result = np.array(cpt_ha_info) @ np.array(HA_COEF_COEF) + HA_INCP
    cpt_ha_info['ha_pred'] = cpt_ha_result
    cpt_ha_info['HA'] = cpt_ha_result
    cpt_ha_info['HA'][cpt_ha_info['HA'] >= HA_TH] = 2
    cpt_ha_info['HA'][cpt_ha_info['HA'] < HA_TH] = 0
    total_info = cpt_col_info.join(cpt_ha_info)
    total_info = total_info[['col_pred', 'ha_pred', 'COL', 'HA']]
    total_info['COMB'] = total_info['COL'] + total_info['HA']
    total_info['COL'][total_info['COL'] == 1] = 'COL_H'
    total_info['COL'][total_info['COL'] == 0] = 'COL_L'
    total_info['HA'][total_info['HA'] == 2] = 'HA_H'
    total_info['HA'][total_info['HA'] == 0] = 'HA_L'
    total_info['COMB'][total_info['COMB'] == 0] = 'L'
    total_info['COMB'][total_info['COMB'] == 1] = 'COL'
    total_info['COMB'][total_info['COMB'] == 2] = 'HA'
    total_info['COMB'][total_info['COMB'] == 3] = 'D'
    total_info.to_csv(CFG.lassoresultpath + 'LASSO_CPT_COMB.csv')
