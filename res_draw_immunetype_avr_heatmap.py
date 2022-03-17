#!/usr/bin/env python
# -*- coding: utf-8 -*-
# @Author : Ji Hongchen
# @Email : jhca.xyt@163.com
# @Last version : 2022-01-17 00:11
# @Filename : res_draw_immunetype_avr_heatmap.py

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns

immune_type = pd.read_csv(
    '~/Documents/MANU/BI/PANC_ECM/DATA/Immue_type_pathway.csv')

immune_gsva = pd.read_csv(
    '~/Documents/MANU/BI/PANC_ECM/DATA/Immune_type_gsva_result.csv', index_col=[0])
for i in immune_gsva.index:
    gsva_d_m = immune_gsva.loc[i].mean()
    gsva_d_s = immune_gsva.loc[i].std()
    immune_gsva.loc[i] = (immune_gsva.loc[i] - gsva_d_m) / gsva_d_s
clical_data = pd.read_csv('~/Documents/MANU/BI/PANC_ECM/Manuscript/Data/Clinica_TCGA_TYPE2_COEF_2.csv',
                          index_col=[0])
immune_gsva = immune_gsva.join(clical_data[['lasso_result']])
print(immune_gsva)

for c in immune_type.columns:
    pathway_list = list(set(immune_type[c]))
    if 'Blank' in pathway_list:
        pathway_list.remove('Blank')
    index_list = ['STROMA_HIGH', 'STROMA_LOW']
    drawframe = pd.DataFrame(columns=pathway_list, index=index_list)
    for p in pathway_list:
        for i in index_list:
            drawframe[p][i] = immune_gsva[p][immune_gsva['lasso_result'] == i].mean()
    drawframe = drawframe.astype('float64')
    figure = plt.figure(figsize=(7, 7), dpi=300)
    sns.heatmap(drawframe,
                cbar=True,
                vmin=-0.5,
                vmax=0.5,
                cmap='bwr',
                xticklabels=pathway_list,
                yticklabels=index_list)
    plt.xticks(size=16)
    plt.yticks(size=16)
    plt.gcf().subplots_adjust(left=0.08, right=0.92, bottom=0.75)
    plt.savefig('/Users/freud/Documents/MANU/BI/PANC_ECM/Manuscript/Figure/SubFigure/' +
                'Immune_type_' + c + '_heatmap.tiff', transparent=True)
    plt.close()
