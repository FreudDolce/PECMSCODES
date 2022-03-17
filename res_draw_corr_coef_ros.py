#!/usr/bin/env python
# -*- coding: utf-8 -*-
# @Author : Ji Hongchen
# @Email : jhca.xyt@163.com
# @Last version : 2022-01-12 18:26
# @Filename : res_draw_corr.py

import pandas as pd
import numpy as np
import os
import matplotlib.pyplot as plt
import seaborn as sns

roslist = list(pd.read_csv(
    '~/Documents/MANU/BI/PANC_ECM/Manuscript/Data/ROS_genes.csv', index_col=[0])['ROS'])

coeflist = list(pd.read_csv(
    '~/Documents/MANU/BI/PANC_ECM/Manuscript/Data/TYPE2_lasso_coef2.csv',
    index_col=[0]).index)
coeflist.remove('(Intercept)')

EXP = pd.read_csv(
    '~/Documents/MANU/BI/PANC_ECM/DATA/SELECTED_MRNA_SYMBOL.csv',
    index_col=[0])

EXP = pd.DataFrame(np.array(EXP).T, index=EXP.columns, columns=EXP.index)
EXP = np.log2(EXP)
corrframe = pd.DataFrame(columns=roslist, index=coeflist)

for r in roslist:
    for c in coeflist:
        ros = np.array(EXP[r])
        coef = np.array(EXP[c])
        corr_list = np.corrcoef(ros, coef)
        corr_v = corr_list[0][1]
        corrframe[r][c] = corr_v
corrframe.to_csv('~/Documents/MANU/BI/PANC_ECM/Manuscript/Data/Corr_ros_coef.csv')

# Draw correlationship map
for c in coeflist:
    for r in roslist:
        figure = plt.figure(figsize=(1, 1), dpi=300)
        sns.jointplot(x=c, y=r, data=EXP, kind='reg',
                      size=2.5, height=1)
        plt.title('R2: ' + str(round(corrframe[r][c], 2)),
                  loc='left', horizontalalignment='left', fontsize=12,
                  verticalalignment='bottom',
                  bbox={'facecolor': 'white', 'edgecolor': 'white', 'alpha': 0.5})
        plt.xlabel('')
        plt.ylabel('')
        plt.tight_layout()
        plt.savefig('/Users/freud/Documents/MANU/BI/PANC_ECM/Manuscript/Figure/SubFigure/' +
                    'Coor_' + c + '_' + r + '.png', transparent=True)
        plt.close()
