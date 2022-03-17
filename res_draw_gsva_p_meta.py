#!/usr/bin/env python
# -*- coding: utf-8 -*-
# @Author : Ji Hongchen
# @Email : jhca.xyt@163.com
# @Last version : 2022-01-13 22:48
# @Filename : res_draw_gsva_score_p_value.py

import pandas as pd
import numpy as np
from scipy import stats
import matplotlib.pyplot as plt

gsva_data = pd.read_csv(
    '~/Documents/MANU/BI/PANC_ECM/DATA/gsva_result.csv', index_col=[0])
gsva_data = pd.DataFrame(np.array(gsva_data).T,
                         columns=gsva_data.index,
                         index=gsva_data.columns)

gsva_list = list(pd.read_csv('~/Documents/MANU/BI/PANC_ECM/Manuscript/Data/GSVA_meta_pathway.csv',
                             index_col=[0])['Pathway'])

gsva_data = gsva_data[gsva_list]

class_list = pd.read_csv('~/Documents/MANU/BI/PANC_ECM/Manuscript/Data/Clinica_TCGA_TYPE2_COEF_2.csv',
                         index_col=[0])[['lasso_result']]

gsva_data = gsva_data.join(class_list)

p_stat = pd.DataFrame(columns=['pvalue'])

for pathway in gsva_list:
    exp_high = gsva_data[pathway][gsva_data['lasso_result'] == 'STROMA_HIGH']
    exp_low = gsva_data[pathway][gsva_data['lasso_result'] == 'STROMA_LOW']
    p = stats.ttest_ind(exp_high, exp_low)
    if exp_high.mean() >= exp_low.mean():
        p_stat.loc[pathway] = abs(np.log10(p.pvalue))
    else:
        p_stat.loc[pathway] = -abs(np.log10(p.pvalue))
    p_stat.sort_values(by='pvalue', ascending=True, inplace=True)
print(p_stat)
figure = plt.figure(figsize=(7, 8), dpi=300)
ax = plt.axes()
plt.barh(y=p_stat.index,
         width=p_stat['pvalue'],
         color='red',
         height=0.4)
plt.tight_layout()
ax.spines['top'].set_visible(False)
ax.spines['right'].set_visible(False)
plt.savefig('/Users/freud/Documents/MANU/BI/PANC_ECM/Manuscript/Figure/SubFigure/Meta_pathway_gsva.tiff')
