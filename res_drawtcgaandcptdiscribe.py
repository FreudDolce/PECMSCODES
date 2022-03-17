#!/usr/bin/env python
# -*- coding: utf-8 -*-
# @Author : Ji Hongchen
# @Email : jhca.xyt@163.com
# @Last version : 2022-01-05 21:03
# @Filename : res_drawtcgaandcptdiscribe.py

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import cfg

item = 'TYPE2'
COEF = '2'


def DrawBoxMap(class_by, value, dataframe, lw, draw_order=[]):
    figure = plt.figure(figsize=(2.6, 8), dpi=300)
    sns.boxplot(x=class_by, y=value,
                data=dataframe, color='white',
                linewidth=lw,
                width=0.7,
                order=draw_order)
    plt.xticks(rotation='vertical', size=20)
    plt.yticks(size=20)
    return plt

cpt_file = pd.read_csv('~/Documents/MANU/BI/PANC_ECM/Manuscript/Data/Clinica_CPT_' + item + '_COEF_' + COEF + '.csv', index_col=[0])
tcga_file = pd.read_csv('~/Documents/MANU/BI/PANC_ECM/Manuscript/Data/Clinica_TCGA_' + item + '_COEF_' + COEF + '.csv', index_col=[0])

drawframe = pd.DataFrame(columns=['class', 'pred'])

for i in cpt_file.index:
    drawframe.loc[i] = list(cpt_file[['project_id', 'lasso_pred']].loc[i])
for j in tcga_file.index:
    drawframe.loc[j] = list(tcga_file[['project_id', 'lasso_pred']].loc[j])

print(drawframe)
drawframe.to_csv('~/Documents/MANU/BI/PANC_ECM/Manuscript/Data/cpt_tcga_pred' + item + '_COEF_' + COEF + '.csv')

plt = DrawBoxMap(class_by='class',
                 value='pred',
                 dataframe=drawframe,
                 lw=2.5,
                 draw_order=['CPTAC-3', 'TCGA-PAAD'])
sns.swarmplot(x='class', y='pred',
              data=drawframe,
              palette='deep',
              order=['CPTAC-3', 'TCGA-PAAD'],
              size=7)
plt.tight_layout()
plt.savefig('/Users/freud/Documents/MANU/BI/PANC_ECM/Manuscript/Figure/SubFigure/CPT_TCGA_double_pred_scatter_'
            + item + '_COEF_' + COEF + '.tiff')
plt.close()
