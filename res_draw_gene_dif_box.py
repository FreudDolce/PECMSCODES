#!/usr/bin/env python
# -*- coding: utf-8 -*-
# @Author : Ji Hongchen
# @Email : jhca.xyt@163.com
# @Last version : 2022-01-12 23:09
# @Filename : res_draw_gene_dif_box.py

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from scipy import stats

Drawlist = list(pd.read_csv('~/Documents/MANU/BI/PANC_ECM/Manuscript/Data/Corr_ros_coef.csv',
                            index_col=[0]).columns)

EXP = pd.read_csv('~/Documents/MANU/BI/PANC_ECM/DATA/SELECTED_MRNA_SYMBOL.csv',
                  index_col=[0])
EXP = pd.DataFrame(np.array(EXP).T,
                   columns=EXP.index,
                   index=EXP.columns)
EXP = np.log2(EXP)
exp_m = np.array(EXP).mean()
exp_s = np.array(EXP).std()


CPT_EXP = pd.read_csv('~/Documents/MANU/BI/PANC_ECM/DATA/CPTAC_MRNA_SYMBOL.csv',
                      index_col=[0])
CPT_EXP = np.log2(CPT_EXP)
CPT_EXP.replace(-np.inf, np.nan, inplace=True)
for cpt_e in CPT_EXP:
    CPT_EXP[cpt_e].fillna(CPT_EXP[cpt_e].mean(), inplace=True)
cpt_m = np.array(CPT_EXP).mean()
cpt_s = np.array(CPT_EXP).std()

#CPT_EXP = (CPT_EXP - cpt_m) / cpt_s
#CPT_EXP = CPT_EXP * exp_s + exp_m

clic = pd.read_csv('~/Documents/MANU/BI/PANC_ECM/Manuscript/Data/Clinica_COMB_TYPE2_COEF_2.csv',
                   index_col=0)[['project_id', 'lasso_result']]
EXP = pd.concat((EXP, CPT_EXP), axis=0, join='inner', sort=True)
for i in EXP.index:
    exp_i_m = EXP.loc[i].mean()
    exp_i_s = EXP.loc[i].std()
    EXP.loc[i] = (EXP.loc[i] - exp_i_m) / exp_i_s
print(EXP)
EXP = EXP[Drawlist]
EXP = EXP.join(clic)


for gene in Drawlist:
    gene_in_high = np.array(EXP[gene][EXP['lasso_result'] == 'STROMA_HIGH'])
    gene_in_low = np.array(EXP[gene][EXP['lasso_result'] == 'STROMA_LOW'])
    p = stats.ttest_ind(gene_in_high, gene_in_low)
    p = stats.mannwhitneyu(gene_in_high, gene_in_low)
    print(p.pvalue)

    figure = plt.figure(figsize=(2.35, 4), dpi=300)
    sns.swarmplot(x='lasso_result',
                  y=gene,
                  hue='project_id',
                  palette='deep',
                  alpha=0.8,
                  data=EXP,
                  size=3.5,
                  order=['STROMA_HIGH', 'STROMA_LOW']
                  )
    sns.boxplot(x='lasso_result',
                y=gene,
                data=EXP,
                #palette='deep',
                color='white',
                linewidth=3,
                width=0.5,
                saturation=1,
                order=['STROMA_HIGH', 'STROMA_LOW']
                )
    plt.title('P = ' + str(round(p.pvalue, 6)))
    bwith = 2.5 #边框宽度设置为2
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
    plt.gcf().subplots_adjust(left=0.4, right=0.92)
    # plt.tight_layout()
    plt.savefig(
        '/Users/freud/Documents/MANU/BI/PANC_ECM/Manuscript/Figure/SubFigure/Comb_ROS_' +
        gene + '_box.tiff', transparent=True)
