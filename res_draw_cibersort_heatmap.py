#!/usr/bin/env python
# -*- coding: utf-8 -*-
# @Author : Ji Hongchen
# @Email : jhca.xyt@163.com
# @Last version : 2022-01-18 21:38
# @Filename : res_draw_cibersort_heatmap.py

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from scipy import stats

Ciberdata = pd.read_csv('~/Documents/MANU/BI/PANC_ECM/DATA/CIBERSORT_RESULT_for_plot.csv',
                        index_col=[0])
delitem = ['P-value', 'Correlation', 'RMSE']
for i in delitem:
    del Ciberdata[i]
Drawlist = list(Ciberdata.columns)
clinidata = pd.read_csv('~/Documents/MANU/BI/PANC_ECM/Manuscript/Data/Clinica_TCGA_TYPE2_COEF_2.csv',
                        index_col=[0])

Ciberdata = Ciberdata.join(clinidata)

Drawframe = pd.DataFrame(columns=['p'])
ciber_high = Ciberdata[Drawlist][Ciberdata['lasso_result'] == 'STROMA_HIGH']
ciber_low = Ciberdata[Drawlist][Ciberdata['lasso_result'] == 'STROMA_LOW']

for cell in Drawlist:
    cell_high = ciber_high[cell]
    cell_low = ciber_low[cell]
    p = stats.mannwhitneyu(cell_high, cell_low)
    if cell_high.mean() >= cell_low.mean():
        Drawframe.loc[cell] = abs(np.log10(p.pvalue))
    else:
        Drawframe.loc[cell] = -abs(np.log10(p.pvalue))

Drawframe.sort_values(by='p', ascending=True, inplace=True)
print(Drawframe)
figure = plt.figure(figsize=(5.5, 10), dpi=300)
ax = plt.axes()
for i in Drawframe.index:
    if abs(Drawframe['p'][i]) <= abs(np.log10(0.05)):
        plt.barh(y=i,
                 width=Drawframe['p'][i],
                 color='grey',
                 height=0.6)
    elif Drawframe['p'][i] >= np.log10(0.05):
        plt.barh(y=i,
                 width=Drawframe['p'][i],
                 color='r',
                 height=0.6)
    else:
        plt.barh(y=i,
                 width=Drawframe['p'][i],
                 color='b',
                 height=0.6)
plt.axvline(x=0, ls='-', lw=2, c='black')
plt.axvline(x=np.log10(0.05), ls=':', lw=2, c='black')
plt.axvline(x=np.log10(0.0005), ls=':', lw=2, c='black')
plt.axvline(x=-np.log10(0.05), ls=':', lw=2, c='black')
plt.axvline(x=-np.log10(0.0005), ls=':', lw=2, c='black')
bwith = 0 #边框宽度设置为2
ax = plt.gca()#获取边框
ax.spines['bottom'].set_linewidth(bwith)
ax.spines['left'].set_linewidth(bwith)
ax.spines['top'].set_linewidth(bwith)
ax.spines['right'].set_linewidth(bwith)
plt.tick_params(width=2.5, length=5.5)
plt.xlabel(' ')
plt.ylabel(' ')
plt.xticks([])
plt.yticks(fontproperties='Arial', size=25)
plt.legend(bbox_to_anchor=(1, 1.4))
plt.gcf().subplots_adjust(left=0.7, right=0.92)
# plt.tight_layout()
plt.savefig(
    '/Users/freud/Documents/MANU/BI/PANC_ECM/Manuscript/Figure/SubFigure/Ciber_b_bar.tiff', transparent=True)
