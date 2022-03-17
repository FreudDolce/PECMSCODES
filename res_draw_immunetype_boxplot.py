#!/usr/bin/env python
# -*- coding: utf-8 -*-
# @Author : Ji Hongchen
# @Email : jhca.xyt@163.com
# @Last version : 2022-01-16 23:00
# @Filename : res_draw_immune_boxplot.py

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from scipy import stats

EXP = pd.read_csv(
    '~/Documents/MANU/BI/PANC_ECM/DATA/Immune_type_gsva_result.csv', index_col=[0])
Drawlist = list(EXP.columns)
for d in EXP.index:
    EXP_d_m = EXP.loc[d].mean()
    EXP_d_s = EXP.loc[d].std()
    EXP.loc[d] = (EXP.loc[d] - EXP_d_m) / EXP_d_s
print(EXP)


class_list = pd.read_csv('~/Documents/MANU/BI/PANC_ECM/Manuscript/Data/Clinica_TCGA_TYPE2_COEF_2.csv',
                         index_col=[0])[['lasso_result']]

EXP = EXP.join(class_list)

for gene in Drawlist:
    gene_in_high = np.array(EXP[gene][EXP['lasso_result'] == 'STROMA_HIGH'])
    gene_in_low = np.array(EXP[gene][EXP['lasso_result'] == 'STROMA_LOW'])
    p = stats.ttest_ind(gene_in_high, gene_in_low)
    p = stats.mannwhitneyu(gene_in_high, gene_in_low)
    print(p.pvalue)

    figure = plt.figure(figsize=(2.35, 4), dpi=300)
    sns.boxplot(x='lasso_result',
                y=gene,
                data=EXP,
                palette='deep',
                linewidth=2.5,
                width=0.5,
                saturation=1,
                order=['STROMA_HIGH', 'STROMA_LOW']
                )
    plt.title('P = ' + str(round(p.pvalue, 6)))
    bwith = 2.5  # 边框宽度设置为2
    ax = plt.gca()  # 获取边框
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
    plt.gcf().subplots_adjust(left=0.4, right=0.92)
    plt.savefig('/Users/freud/Documents/MANU/BI/PANC_ECM/Manuscript/Figure/SubFigure/Immune_gsva_' + gene + '.tiff', transparent=True)
