#!/usr/bin/env python
# -*- coding: utf-8 -*-
# @Author : Ji Hongchen
# @Email : jhca.xyt@163.com
# @Last version : 2022-01-19 23:05
# @Filename : res_draw_gene_heatmap.py

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


def _Draw_Heat_Map(dataframe, v_min, v_max):
    figure = plt.figure(figsize=(len(dataframe) + 1, 5), dpi=300)
    sns.heatmap(dataframe,
                cbar=True,
                vmin=v_min,
                vmax=v_max,
                xticklabels=dataframe.columns,
                yticklabels=dataframe.index,
                cmap='bwr'
                )
    ax = plt.gca()
    bwith = 3  # 边框宽度设置为2
    ax = plt.gca()  # 获取边框
    ax.spines['bottom'].set_linewidth(bwith)
    ax.spines['left'].set_linewidth(bwith)
    ax.spines['top'].set_linewidth(bwith)
    ax.spines['right'].set_linewidth(bwith)
    plt.tick_params(width=2.5, length=5.5)
    plt.yticks(fontproperties='Arial', size=25)
    plt.xticks(fontproperties='Arial', size=25, rotation=90)
    plt.legend(bbox_to_anchor=(1, 1.4))
    plt.gcf().subplots_adjust(left=0.4, right=0.92)
    return plt


EXP = pd.DataFrame(np.array(EXP).T, columns=EXP.index, index=EXP.columns)

genelist = list(pd.read_csv('~/Documents/MANU/BI/PANC_ECM/Manuscript/Data/TYPE2_lasso_coef2.csv',
                            index_col=[0]).index)
genelist.remove('(Intercept)')
EXP = EXP[genelist]

# Survival heatmap
single_Survframe = pd.read_csv('~/Documents/MANU/BI/PANC_ECM/Manuscript/Data/SingleCoxSig.csv',
                               index_col=[0])[['coef']]
multip_Survframe = pd.read_csv('~/Documents/MANU/BI/PANC_ECM/Manuscript/Data/Part_MutiCox.csv',
                               index_col=[0])[['coef']]
multip_Survframe.columns = ['sig_coef']
Survframe = single_Survframe.join(multip_Survframe)
plt = _Draw_Heat_Map(Survframe,
                     v_min=0,
                     v_max=0.3)
plt.savefig('/Users/freud/Documents/MANU/BI/PANC_ECM/Manuscript/Figure/SubFigure/gene_suvr_cox.tiff',
            transparent=True)
plt.close()

# pathway heatmap
pathway_data = pd.read_csv('~/Documents/MANU/BI/PANC_ECM/DATA/Immune_type_gsva_result.csv',
                           index_col=[0])
pathwaylist = ['', '', '', '', '']
