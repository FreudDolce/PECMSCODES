#!/usr/bin/env python
# -*- coding: utf-8 -*-
# @Author : Ji Hongchen
# @Email : jhca.xyt@163.com
# @Last version : 2022-01-12 22:38
# @Filename : res_draw_coef_ros_heatmap.py

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns

coefl = {'PDCD1': 5133,
         'CD274': 29126,
         'PDCD1LG2': 80380,
         'CTLA4': 1493,
         'CD80': 941,
         'CD86': 942,
         'HAVCR2': 84868,
         'TIGIT': 201633,
         'TNFRSF9': 3604,
         'CD27': 939,
         'TOX': 9760,
         'ENTPD1': 953}
coeflist = {}
for c in coefl:
    coeflist[coefl[c]] = c
print(coeflist)

Heatdata = pd.read_csv('~/Documents/MANU/BI/PANC_ECM/Manuscript/Data/Corr_immune_coef.csv',
                       index_col=[0])

X = []
print(Heatdata.columns)
for x in Heatdata.columns:
    X.append(coeflist[int(x)])
Y = list(Heatdata.index)

figure = plt.figure(figsize=(4.2, 3.1), dpi=300)
plt.style.use('ggplot')
sns.set_style('whitegrid')

sns.heatmap(Heatdata,
            cbar=True,
            vmin=-0.6,
            vmax=0.6,
            cmap='bwr',
            xticklabels=X,
            yticklabels=Y
            )
plt.xticks(size=16.5)
plt.yticks(size=16.5)
plt.tight_layout()
plt.savefig(
    '/Users/freud/Documents/MANU/BI/PANC_ECM/Manuscript/Figure/SubFigure/Coef_immune_corr.tiff')
