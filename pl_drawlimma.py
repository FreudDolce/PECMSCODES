#!/usr/bin/env python
# -*- coding: utf-8 -*-
# @Author : Ji Hongchen
# @Email : jhca.xyt@163.com
# @Last version : 2021-12-26 11:53
# @Filename : pl_drawlimma.py

import pandas as pd
import numpy as np
import os
import matplotlib.pyplot as plt
import cfg

CFG = cfg.cfg()

P_LEVEL = 5
FDR_LEVEL = 5
FC_LEVEL = 1.2
ITEM = 'HA'
PLOT_SIZE = 3
DRAW_FILE = CFG.lassoresultpath + 'SCORE/' + \
    ITEM + '_H-' + ITEM + '_L_exp_dif.csv'

df = pd.read_csv(DRAW_FILE)
df['lp'] = -np.log10(df['p'])
df['lfdr'] = -np.log10(df['fdr'])
rankdf = df[df['lp'] > P_LEVEL]
rankdf = rankdf[rankdf['log2fc'].abs() > FC_LEVEL]
rankdf = rankdf.sort_values(by='lp', ascending=False)[:5]
fdrdf = df[df['lfdr'] > FDR_LEVEL]
fdrdf = fdrdf[fdrdf['log2fc'].abs() > FC_LEVEL]
fdrdf.to_csv(CFG.lassoresultpath + 'SCORE/' + ITEM + '_H-' + ITEM + '_L_sig_dif.csv')
print(rankdf)

df_s = df[df['lp'] > P_LEVEL]
df_r = df_s[df_s['log2fc'] > FC_LEVEL]
df_b = df_s[df_s['log2fc'] < -FC_LEVEL]
X1 = df['log2fc']
Y1 = df['lp']
X2 = df_r['log2fc']
Y2 = df_r['lp']
X3 = df_b['log2fc']
Y3 = df_b['lp']


plt.figure(figsize=(6, 6), dpi=300)

plt.scatter(X1, Y1, s=PLOT_SIZE, c='silver')
plt.scatter(X2, Y2, s=PLOT_SIZE, c='orangered')
plt.scatter(X3, Y3, s=PLOT_SIZE, c='dodgerblue')
for point in rankdf.index:
    if rankdf['log2fc'].loc[point] > 0:
        plt.text(
            rankdf['log2fc'].loc[point],
            rankdf['lp'].loc[point],
            rankdf['X'].loc[point],
            fontsize=7,
            rotation=-45,
            color='orangered',
            verticalalignment='top',
            horizontalalignment='left')
    if rankdf['log2fc'].loc[point] < 0:
        plt.text(
            rankdf['log2fc'].loc[point],
            rankdf['lp'].loc[point],
            rankdf['X'].loc[point],
            fontsize=7,
            rotation=-45,
            color='dodgerblue',
            verticalalignment='top',
            horizontalalignment='left')
plt.axhline(y=P_LEVEL, ls=':', c='gray', lw=1)
plt.axvline(x=FC_LEVEL, ls=':', c='gray', lw=1)
plt.axvline(x=-FC_LEVEL, ls=':', c='gray', lw=1)
plt.xlim(-3, 3)
plt.tight_layout()

plt.savefig(CFG.lassoresultpath + ITEM + '_LIMMA_PLOT.jpg')
