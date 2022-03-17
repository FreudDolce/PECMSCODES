#!/usr/bin/env python
# -*- coding: utf-8 -*-
# @Author : Ji Hongchen
# @Email : jhca.xyt@163.com
# @Last version : 2022-01-21 00:11
# @Filename : res_draw_gene2gene_corr.py

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns

EXP = pd.read_csv('~/Documents/MANU/BI/PANC_ECM/DATA/SELECTED_MRNA_SYMBOL.csv',
                  index_col=[0])
EXP = np.log2(EXP)
for c in EXP.columns:
    exp_c_m = EXP[c].mean()
    exp_c_s = EXP[c].std()
    EXP[c] = (EXP[c] - exp_c_m) / exp_c_s

EXP = pd.DataFrame(np.array(EXP).T, index=EXP.columns, columns=EXP.index)

genelist = list(pd.read_csv('~/Documents/MANU/BI/PANC_ECM/Manuscript/Data/TYPE2_lasso_coef2.csv',
                            index_col=[0]).index)
genelist.remove('(Intercept)')

EXP = EXP[genelist]

Drawframe = pd.DataFrame(columns=genelist, index=genelist)
for g in genelist:
    for k in genelist:
        corr_v = np.corrcoef(EXP[g], EXP[k])[0][1]
        Drawframe[g][k] = corr_v

Drawframe = Drawframe.astype('float')
print(Drawframe)
figure = plt.figure(figsize=(6, 6), dpi=300)
for i in range(len(genelist)):
    for j in range(i, len(genelist)):
        corr = Drawframe[genelist[i]][genelist[j]]
        print(corr)
        plt.scatter(x=[i],
                    y=[j],
                    s=abs(corr * 20) ** 2,
                    color=((100 * corr + 130) / 255,
                           0,
                           1.1 - (100 * corr + 130) / 255)
                    )

bwith = 0  # 边框宽度设置为2
ax = plt.gca()  # 获取边框
ax.spines['bottom'].set_linewidth(bwith)
ax.spines['left'].set_linewidth(bwith)
ax.spines['top'].set_linewidth(bwith)
ax.spines['right'].set_linewidth(bwith)
plt.tick_params(width=2.5, length=5.5)
plt.xlabel(' ')
plt.ylabel(' ')
plt.yticks(np.arange(len(genelist)),
           genelist,
           fontproperties='Arial', size=25)
plt.xticks(np.arange(len(genelist)),
           genelist,
           fontproperties='Arial', size=25)
plt.legend(bbox_to_anchor=(1, 1.4))
plt.gcf().subplots_adjust(left=0.4, right=0.92, bottom=0.4, top=0.92)
# plt.tight_layout()
plt.savefig(
    '/Users/freud/Documents/MANU/BI/PANC_ECM/Manuscript/Figure/SubFigure/Comb_gene2gene_corr.tiff')
