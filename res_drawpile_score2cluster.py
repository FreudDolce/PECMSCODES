#!/usr/bin/env python
# -*- coding: utf-8 -*-
# @Author : Ji Hongchen
# @Email : jhca.xyt@163.com
# @Last version : 2022-01-05 13:59
# @Filename : res_drawpile_score2cluster.py

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import cfg

db = 'TCGA'
item = 'TYPE2'
COEF = '2'
CFG = cfg.cfg()
colors = CFG.colors


def DrawPileMap(dataframe, class_col, draw_col, draw_order='nan'):
    hue_class = list(set(dataframe[class_col]))
    hue_class.sort()
    draw_class = list(set(dataframe[draw_col]))
    bottom_list = [0] * len(hue_class)
    figure = plt.figure(figsize=(len(hue_class)*1.3 + 1, 10), dpi=300)
    if draw_order == 'nan':
        for d in range(len(draw_class)):
            drawlist = []
            for h in hue_class:
                hue = dataframe[dataframe[class_col] == h]
                draw_hue = len(hue[hue[draw_col] == draw_class[d]]) / len(hue)
                drawlist.append(draw_hue)
            plt.bar(hue_class, drawlist, bottom=bottom_list,
                    label=draw_class[d], color=colors[d])
            for i in range(len(drawlist)):
                bottom_list[i] += drawlist[i]
    else:
        print('Draw with order!')
        for d in range(len(draw_order)):
            drawlist = []
            for h in hue_class:
                hue = dataframe[dataframe[class_col] == h]
                draw_hue = len(hue[hue[draw_col] == draw_order[d]]) / len(hue)
                drawlist.append(draw_hue)
            plt.bar(hue_class, drawlist, bottom=bottom_list,
                    label=draw_order[d], color=colors[d])
            for i in range(len(drawlist)):
                bottom_list[i] += drawlist[i]
    plt.ylim(-0.02, 1.02)
    plt.xticks(range(len(hue_class)), hue_class, rotation='vertical', size=20)
    plt.yticks([0, 0.25, 0.5, 0.75, 1],
               ['0%', '25%', '50%', '75%', '100%'],
               size=20)
    plt.legend(loc=[0, 1.1])
    return plt


def DrawBoxMap(class_by, value, dataframe, lw, draw_order=[]):
    figure = plt.figure(figsize=(4.5, 8), dpi=300)
    sns.boxplot(x=class_by, y=value,
                data=dataframe, color='white',
                linewidth=lw,
                width=0.7,
                order=draw_order)
    plt.xticks(rotation='vertical', size=20)
    plt.yticks(size=20)
    return plt


cluster_order = ['Cluster_1', 'Cluster_2', 'Cluster_3', 'Cluster_4']

score2cluster = pd.read_csv('~/Documents/MANU/BI/PANC_ECM/Manuscript/Data/' +
                            'Clinica_' + db + '_' + item + '_COEF_' + COEF + '.csv',
                            index_col=[0])

plt = DrawPileMap(score2cluster, 'lasso_result',
                  'cluster_result', draw_order=cluster_order)
plt.tight_layout()
plt.savefig('/Users/freud/Documents/MANU/BI/PANC_ECM/Manuscript/Figure/SubFigure/cluster_score_pile_' +
            db + '_' + item + '_COEF_' + COEF + '.tiff')

plt.close()

plt = DrawBoxMap(class_by='cluster_result',
                 value='lasso_pred',
                 dataframe=score2cluster,
                 lw=2.5,
                 draw_order=cluster_order)
sns.swarmplot(x='cluster_result', y='lasso_pred',
              data=score2cluster, palette='deep',
              order=cluster_order,
              size=7)
plt.axhline(0.434, ls=':', c='black', lw=3)
plt.ylim(-0.2, 0.8)
plt.tight_layout()
plt.savefig('/Users/freud/Documents/MANU/BI/PANC_ECM/Manuscript/Figure/SubFigure/score_cluster_box_' +
            db + '_' + item + '_COEF_' + COEF + '.tiff')
plt.close()
# plt = DrawStripplot(dataframe=score2cluster,
#                    class_by='cluster_result',
#                    value='lasso_pred',
#                    draw_order=cluster_order)
# plt.show()
