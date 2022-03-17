#!/usr/bin/env python
# -*- coding: utf-8 -*-
# @Author : Ji Hongchen
# @Email : jhca.xyt@163.com
# @Last version : 2022-01-11 19:49
# @Filename : res_clinical_feature.py

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import cfg
CFG = cfg.cfg()
colors = CFG.colors

REP = {'ajcc_pathologic_stage': {'I': 0, 'IA': 0, 'IB': 0, 'IIA': 0, 'IIB': 0,
                                 'III': 1, 'IIIB': 1, 'IV': 1}}

cpt_clic = pd.read_csv('~/Documents/MANU/BI/PANC_ECM/Manuscript/Data/Clinica_CPT_TYPE2_COEF_2.csv',
                       index_col=[0])
for item in REP['ajcc_pathologic_stage']:
    cpt_clic['ajcc_pathologic_stage'][cpt_clic['ajcc_pathologic_stage'] == item] = \
        REP['ajcc_pathologic_stage'][item]
cpt_clic['days_to_birth'][cpt_clic['days_to_birth'] < 65] = 0
cpt_clic['days_to_birth'][cpt_clic['days_to_birth'] >= 65] = 1

tcga_clic = pd.read_csv('~/Documents/MANU/BI/PANC_ECM/Manuscript/Data/Clinica_TCGA_TYPE2_COEF_2.csv',
                        index_col=[0])
for item in REP['ajcc_pathologic_stage']:
    tcga_clic['ajcc_pathologic_stage'][tcga_clic['ajcc_pathologic_stage'] == item] = \
        REP['ajcc_pathologic_stage'][item]
tcga_clic['days_to_birth'][tcga_clic['days_to_birth'] < 65] = 0
tcga_clic['days_to_birth'][tcga_clic['days_to_birth'] >= 65] = 1

cpt_high = cpt_clic[cpt_clic['lasso_result'] == 'STROMA_HIGH']
print(cpt_high)
cpt_low = cpt_clic[cpt_clic['lasso_result'] == 'STROMA_LOW']
tcga_high = tcga_clic[tcga_clic['lasso_result'] == 'STROMA-HIGH']
tcga_low = tcga_clic[tcga_clic['lasso_result'] == 'STROMA-LOW']

cols = ['tcga_high', 'tcga_low', 'cpt_high', 'cpt_low']
ins = np.array(['tissue_or_organ_of_origin', 'ajcc_pathologic_stage',
    'days_to_birth', 'alcohol_history', 'gender', 'cigarettes_per_day'])[: : -1]

draw_frame = pd.DataFrame(columns=cols, index=ins)

for i in ins:
    draw_frame['tcga_high'].loc[i] = len(
        tcga_high[tcga_high[i] == 1]) / len(tcga_high)
    draw_frame['tcga_low'].loc[i] = len(
        tcga_low[tcga_low[i] == 1]) / len(tcga_low)
    draw_frame['cpt_high'].loc[i] = len(
        cpt_high[cpt_high[i] == 1]) / len(cpt_high)
    draw_frame['cpt_low'].loc[i] = len(cpt_low[cpt_low[i] == 1]) / len(cpt_low)

figure = plt.figure(figsize=(5, 5.3), dpi=300)
ax = plt.axes()
ax.spines['top'].set_visible(False)
ax.spines['bottom'].set_visible(False)
ax.spines['right'].set_visible(False)
ax.spines['left'].set_visible(False)
y_tcga = np.arange(1, len(draw_frame) + 1) - 0.5
y_cpt = np.arange(1, len(draw_frame) + 1) - 0.1
base_pos = np.ones(len(draw_frame))
base_neg = -1 * np.ones(len(draw_frame))
plt.barh(y=y_tcga,
         width=base_pos,
         height=0.4,
         color=colors[0],
         alpha=0.3)
plt.barh(y=y_cpt,
         width=base_pos,
         height=0.4,
         color=colors[1],
         alpha=0.3)
plt.barh(y=y_tcga,
         width=base_neg,
         height=0.4,
         color=colors[0],
         alpha=0.3)
plt.barh(y=y_cpt,
         width=base_neg,
         height=0.4,
         color=colors[1],
         alpha=0.3)
plt.barh(y=y_tcga,
         width=draw_frame['tcga_high'],
         height=0.4,
         color=colors[0])
plt.barh(y=y_cpt,
         width=draw_frame['cpt_high'],
         height=0.4,
         color=colors[1])
plt.barh(y=y_tcga,
         width=-1 * draw_frame['tcga_low'],
         height=0.4,
         color=colors[0])
plt.barh(y=y_cpt,
         width=-1 * draw_frame['cpt_low'],
         height=0.4,
         color=colors[1])
plt.axvline(x=0, ls=':', lw=2, c='black')
plt.axvline(x=0.25, ls=':', lw=1, c='black')
plt.axvline(x=0.5, ls=':', lw=1, c='black')
plt.axvline(x=0.75, ls=':', lw=1, c='black')
plt.axvline(x=-0.25, ls=':', lw=1, c='black')
plt.axvline(x=-0.5, ls=':', lw=1, c='black')
plt.axvline(x=-0.75, ls=':', lw=1, c='black')
plt.axhline(y=-1, ls='-', lw=2, c='black')
plt.yticks(np.arange(len(draw_frame)) + 0.5, np.array(draw_frame.index))
plt.xticks([-1, -0.5, 0, 0.5, 1], ['100%', '50%', '0%', '50%', '100%'])
plt.ylim(-1.1, 7.1)
plt.tight_layout()
plt.savefig(
    '/Users/freud/Documents/MANU/BI/PANC_ECM/Manuscript/Figure/SubFigure/Clinical_feature.tiff',
    transparent=True)
