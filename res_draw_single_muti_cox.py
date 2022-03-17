#!/usr/bin/env python
# -*- coding: utf-8 -*-
# @Author : Ji Hongchen
# @Email : jhca.xyt@163.com
# @Last version : 2022-01-11 18:36
# @Filename : res_draw_single_muti_cox.py

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

sigle_frame = pd.read_csv('~/Documents/MANU/BI/PANC_ECM/Manuscript/Data/Part_MutiCox.csv',
                          index_col=[0])
sigle_frame.sort_values(by=['coef', 'lower'], ascending=False, inplace=True)
print(sigle_frame)

MUL_COX = pd.read_csv('/Users/freud/Documents/MANU/BI/PANC_ECM/Manuscript/Data/Part_MutiCox.csv',
                      index_col=[0])

multi_frame = pd.DataFrame(
    columns=sigle_frame.columns, index=sigle_frame.index)

for i in MUL_COX.columns:
    for j in MUL_COX.index:
        multi_frame[i].loc[j] = MUL_COX[i][j]

multi_frame.fillna(-3, inplace=True)

y = np.arange(len(sigle_frame))[:: -1] + 0.7

sigle_frame['d_l'] = sigle_frame['coef'] - sigle_frame['lower']
sigle_frame['d_h'] = sigle_frame['upper'] - sigle_frame['coef']
multi_frame['d_l'] = multi_frame['coef'] - multi_frame['lower']
multi_frame['d_h'] = multi_frame['upper'] - multi_frame['coef']

# Single cox:
figure = plt.figure(figsize=(5, 5.3), dpi=300)
ax = plt.axes()
ax.spines['top'].set_visible(False)
ax.spines['bottom'].set_visible(False)
ax.spines['right'].set_visible(False)
plt.errorbar(
    np.array(sigle_frame['coef']),
    y,
    fmt='o',
    color='black',
    ecolor='black',
    elinewidth=2,
    capthick=2,
    capsize=4,
    xerr=np.array(sigle_frame[['d_l', 'd_h']]).T)
plt.xlim(-1.2, 1.2)
plt.ylim(-1.1, 7.5)
plt.yticks(np.arange(len(sigle_frame)) + 0.7,
           np.array(list(sigle_frame.index))[:: -1])
plt.axvline(x=0, ls=':', lw=2, c='black')
plt.axhline(y=-1, ls='-', lw=2, c='black')
plt.tight_layout()
plt.savefig('/Users/freud/Documents/MANU/BI/PANC_ECM/Manuscript/Figure/SubFigure/Gene_SigCox.tiff',
            transparent=True)
plt.close()

# Multiple cox:
figure = plt.figure(figsize=(4, 5.3), dpi=300)
ax = plt.axes()
ax.spines['top'].set_visible(False)
ax.spines['bottom'].set_visible(False)
ax.spines['right'].set_visible(False)
ax.spines['left'].set_visible(False)
plt.errorbar(
    np.array(multi_frame['coef']),
    y,
    fmt='o',
    color='black',
    ecolor='black',
    elinewidth=2,
    capthick=2,
    capsize=4,
    xerr=np.array(multi_frame[['d_l', 'd_h']]).T)
plt.xlim(-1.2, 1.2)
plt.ylim(-1.1, 7.5)
plt.yticks(np.arange(len(multi_frame)) + 0.7,
           np.array(list(multi_frame.index))[:: -1])
plt.axvline(x=0, ls=':', lw=2, c='black')
plt.axhline(y=-1, ls='-', lw=2, c='black')
plt.tight_layout()
plt.savefig(
    '/Users/freud/Documents/MANU/BI/PANC_ECM/Manuscript/Figure/SubFigure/Gene_MulCox.tiff',
    transparent=True)
plt.close()
