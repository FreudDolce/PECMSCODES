#!/usr/bin/env python
# -*- coding: utf-8 -*-
# @Author : Ji Hongchen
# @Email : jhca.xyt@163.com
# @Last version : 2022-01-19 12:25
# @Filename : res_draw_drug_sensitive_boxplot.py

import pandas as pd
import numpy as np
import os
import cfg
from scipy import stats
import matplotlib.pyplot as plt
import seaborn as sns

DrugSensitive = pd.read_csv('~/Documents/MANU/BI/PANC_ECM/DATA/DrugSensitive_predict.csv',
                            index_col=[0])
Drawlist = list(DrugSensitive.columns)

clinicl_data = pd.read_csv('~/Documents/MANU/BI/PANC_ECM/Manuscript/Data/Clinica_COMB_TYPE2_COEF_2.csv',
                           index_col=[0])[['lasso_result', 'lasso_pred']]

DrugSensitive = DrugSensitive.join(clinicl_data, how='inner')

for drug in Drawlist:
    drug_high = np.array(
        DrugSensitive[drug][DrugSensitive['lasso_result'] == 'STROMA_HIGH'])
    drug_low = np.array(
        DrugSensitive[drug][DrugSensitive['lasso_result'] == 'STROMA_LOW'])
    #p = stats.mannwhitneyu(drug_high, drug_low)
    p = stats.ttest_ind(drug_high, drug_low)

    figure = plt.figure(figsize=(2.35, 4), dpi=300)
    sns.boxplot(x='lasso_result',
                y=drug,
                data=DrugSensitive,
                palette='deep',
                # color='white',
                linewidth=3,
                width=0.5,
                saturation=1,
                order=['STROMA_LOW', 'STROMA_HIGH']
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
    # plt.tight_layout()
    plt.savefig(
        '/Users/freud/Documents/MANU/BI/PANC_ECM/Manuscript/Figure/SubFigure/Drug_sensitive_' +
        drug + '_box.tiff', transparent=True)
    plt.close()
    figure = plt.figure(figsize=(1.7, 2), dpi=300)
    corr = np.corrcoef(DrugSensitive[drug], DrugSensitive['lasso_pred'])[0][1]
    sns.lmplot(y=drug,
               x='lasso_pred',
               data=DrugSensitive,
               aspect=0.5)
    plt.title('ρ: ' + str(round(corr, 2)),
              loc='center', horizontalalignment='center', fontsize=25,
              verticalalignment='bottom',
              bbox={'facecolor': 'white', 'edgecolor': 'white', 'alpha': 0.5})
    plt.grid()
    bwith = 2.5  # 边框宽度设置为2
    ax = plt.gca()  # 获取边框
    ax.spines['bottom'].set_linewidth(bwith)
    ax.spines['left'].set_linewidth(bwith)
    ax.spines['top'].set_linewidth(bwith)
    ax.spines['right'].set_linewidth(bwith)
    plt.tick_params(width=2.5, length=5.5)
    plt.xticks(fontproperties='Arial', size=25, rotation=90)
    plt.yticks(fontproperties='Arial', size=25)
    #plt.xlim(-5, 5)
    #plt.ylim(-5, 5)
    plt.xlabel(' ')
    plt.ylabel(' ')
    plt.gcf().subplots_adjust(left=0.4, right=0.92, top=0.8, bottom=0.2)
    # plt.tight_layout()
    plt.savefig(
        '/Users/freud/Documents/MANU/BI/PANC_ECM/Manuscript/Figure/SubFigure/Drug_sensitive_' +
        'Coor_' + drug + '.png', transparent=True)
    plt.close()
