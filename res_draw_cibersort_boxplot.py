#!/usr/bin/env python
# -*- coding: utf-8 -*-
# @Author : Ji Hongchen
# @Email : jhca.xyt@163.com
# @Last version : 2022-01-18 20:41
# @Filename : res_draw_cibersort_boxplot.py

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from scipy import stats

EXP = pd.read_csv(
    '~/Documents/MANU/BI/PANC_ECM/DATA/CIBERSORT_RESULT_for_plot.csv', index_col=[0])
delitem = ['P-value', 'Correlation', 'RMSE']
for i in delitem:
    del EXP[i]
#print(EXP)
#EXP = np.log10(EXP)
#EXP.replace(-np.inf, np.nan, inplace=True)
#EXP.fillna(EXP.median(), inplace=True)
#EXP.fillna(0, inplace=True)
Drawlist = list(set(EXP.columns))

clic = pd.read_csv('~/Documents/MANU/BI/PANC_ECM/Manuscript/Data/Clinica_TCGA_TYPE2_COEF_2.csv',
                   index_col=[0])[['project_id', 'lasso_result']]

EXP = EXP.join(clic)


for gene in Drawlist:
    gene_in_high = np.array(EXP[gene][EXP['lasso_result'] == 'STROMA_HIGH'])
    gene_in_low = np.array(EXP[gene][EXP['lasso_result'] == 'STROMA_LOW'])
    p = stats.ttest_ind(gene_in_high, gene_in_low)
    p = stats.mannwhitneyu(gene_in_high, gene_in_low)
    print(p.pvalue)

    figure = plt.figure(figsize=(1.7, 4), dpi=300)
    sns.boxplot(x='lasso_result',
                y=gene,
                data=EXP,
                palette='deep',
                #color='white',
                linewidth=3,
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
    plt.yticks(fontproperties='Arial', size=15)
    plt.legend(bbox_to_anchor=(1, 1.4))
    plt.gcf().subplots_adjust(left=0.6, right=0.92)
    # plt.tight_layout()
    plt.savefig(
        '/Users/freud/Documents/MANU/BI/PANC_ECM/Manuscript/Figure/SubFigure/Cibersort_' +
        gene + '_box.tiff', transparent=True)
