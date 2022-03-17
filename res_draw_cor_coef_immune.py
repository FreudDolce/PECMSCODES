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

print('PDCD1: 5133; CD274: 29126; PDCD1LG2: 80380; CTLA4: 1493; CD80: 941; CD86: 942; HAVCR2: 84868; TIGIT: 201633; TNFRSF9: 3604; CD27: 939; TOX: 9760; ENTPD1: 953')
roslist = [5133, 29126, 80380, 1493, 941,
           942, 84868, 201633, 3604, 939, 9760, 953]
coeflist = {'COL17A1': 1308, 'AREG': 374, 'KLHL32': 114792,
            'CDA': 978, 'POSTN': 10631, 'SLC2A1': 6513, 'FN1': 2335, 'INHBA': 3624}

EXP = pd.read_csv(
    '~/Documents/MANU/BI/PANC_ECM/DATA/SELECTED_MRNA.csv',
    index_col=[0])

EXP = pd.DataFrame(np.array(EXP).T, index=EXP.columns, columns=EXP.index)
EXP = np.log2(EXP)
EXP.replace(-np.inf, np.nan, inplace=True)
EXP.fillna(0, inplace=True)
corrframe = pd.DataFrame(columns=roslist, index=coeflist)

for r in roslist:
    for c in coeflist:
        ros = np.array(EXP[r])
        coef = np.array(EXP[coeflist[c]])
        corr_list = np.corrcoef(ros, coef)
        corr_v = corr_list[0][1]
        corrframe[r][c] = corr_v
corrframe.to_csv(
    '~/Documents/MANU/BI/PANC_ECM/Manuscript/Data/Corr_immune_coef.csv')

# Draw correlationship map
for c in coeflist:
    for r in roslist:
        figure = plt.figure(figsize=(1, 1), dpi=300)
        sns.jointplot(x=coeflist[c], y=r, data=EXP, kind='reg',
                      size=2.5, height=1)
        plt.title('R2: ' + str(round(corrframe[r][c], 2)),
                  loc='left', horizontalalignment='left', fontsize=12,
                  verticalalignment='bottom',
                  bbox={'facecolor': 'white', 'edgecolor': 'white', 'alpha': 0.5})
        plt.xlabel('')
        plt.ylabel('')
        plt.tight_layout()
        plt.savefig('/Users/freud/Documents/MANU/BI/PANC_ECM/Manuscript/Figure/SubFigure/' +
                    'Immune_coor_' + c + '_' + str(r) + '.png', transparent=True)
        plt.close()
