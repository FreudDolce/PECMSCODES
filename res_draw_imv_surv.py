#!/usr/bin/env python
# -*- coding: utf-8 -*-
# @Author : Ji Hongchen
# @Email : jhca.xyt@163.com
# @Last version : 2022-01-29 23:08
# @Filename : res_draw_imv_surv.py

import pandas as pd
import numpy as np
import os
from res_clusterandclassdraw import DrawPileMap
from scipy.stats import chi2_contingency

thresholdlist = 0.434

imv_exp = pd.read_csv('~/Documents/MANU/BI/PANC_ECM/DATA/IMV_count.txt',
                      sep='\t',
                      index_col=[0])
imv_exp = pd.DataFrame(np.array(imv_exp).T,
                       index=imv_exp.columns, columns=imv_exp.index)

paad_exp = pd.read_csv('~/Documents/MANU/BI/PANC_ECM/DATA/SELECTED_MRNA_SYMBOL.csv',
                       index_col=[0])
paad_exp = pd.DataFrame(np.array(paad_exp).T,
                        index=paad_exp.columns, columns=paad_exp.index)

coeffile = pd.read_csv('~/Documents/MANU/BI/PANC_ECM/Manuscript/Data/TYPE2_lasso_coef1.csv',
                       index_col=[0])

COEF = pd.read_csv(
    '~/Documents/MANU/BI/PANC_ECM/Manuscript/Data/TYPE2_lasso_coef2.csv',
    index_col=[0])
COEF.drop(COEF[COEF['s1'] == 0].index, inplace=True)
INCP = COEF['s1'].loc['(Intercept)']
COEF_LIST = list(COEF.index)
COEF_LIST.remove('(Intercept)')
COEF_COEF = np.array(COEF.loc[COEF_LIST]).reshape(-1, 1)

imv_exp = imv_exp[COEF_LIST]
imv_exp = np.log2(imv_exp)
imv_exp.replace(-np.inf, np.nan, inplace=True)
imv_exp.fillna(0, inplace=True)
paad_exp = paad_exp[COEF_LIST]
paad_exp = np.log2(paad_exp)

imv_m = np.array(imv_exp).mean()
imv_s = np.array(imv_exp).std()

paad_m = np.array(paad_exp).mean()
paad_s = np.array(paad_exp).std()

imv_exp = (imv_exp - imv_m) / imv_s
imv_exp = imv_exp * paad_s + paad_m

imv_result = np.array(imv_exp) @ np.array(COEF_COEF).reshape(-1, 1) + INCP
imv_exp['lasso_pred'] = imv_result
imv_result = imv_exp[['lasso_pred']]

imv_clic = pd.read_csv('~/Documents/MANU/BI/PANC_ECM/DATA/IMV-clinical.txt', sep='\t',
                       index_col=[0])
imv_clic = imv_clic.join(imv_result)
imv_clic['lasso_result'] = 'a'
imv_clic['lasso_result'][imv_clic['lasso_pred']
                         >= thresholdlist] = 'STROMA_HIGH'
imv_clic['lasso_result'][imv_clic['lasso_pred'] < thresholdlist] = 'STROMA_LOW'
# imv_clic.to_csv('~/Documents/MANU/BI/PANC_ECM/Manuscript/Data/IMV_clinical.csv')

rsp_clic = imv_clic.drop(
    imv_clic[imv_clic['Best Confirmed Overall Response'] == 'NE'].index)

rsp_high = rsp_clic[rsp_clic['lasso_result'] == 'STROMA_HIGH']
rsp_low = rsp_clic[rsp_clic['lasso_result'] == 'STROMA_LOW']

DRAW_ORDER = ['CR', 'PR', 'SD', 'PD']
drawframe = pd.DataFrame(columns=['STROMA_HIGH', 'STROMA_LOW'],
                         index=DRAW_ORDER)

for d in DRAW_ORDER:
    drawframe['STROMA_HIGH'].loc[d] = \
            len(rsp_high[rsp_high['Best Confirmed Overall Response'] == d])
    drawframe['STROMA_LOW'].loc[d] = \
            len(rsp_low[rsp_low['Best Confirmed Overall Response'] == d])

drawframe.to_csv(
        '~/Documents/MANU/BI/PANC_ECM/Manuscript/Data/IMV_response.csv')
kdf, p, dof, expctd = chi2_contingency(drawframe)
print(p)


plt = DrawPileMap(rsp_clic,
                  'lasso_result',
                  'Best Confirmed Overall Response',
                  draw_order=DRAW_ORDER)
plt.legend(bbox_to_anchor=(0, 1.2))
plt.yticks(fontproperties='Arial', size=25)
plt.gcf().subplots_adjust(left=0.4, right=0.92, bottom=0.3)
plt.savefig(
    '/Users/freud/Documents/MANU/BI/PANC_ECM/Manuscript/Figure/SubFigure/IMV_response.tiff',
    transparent=True)
