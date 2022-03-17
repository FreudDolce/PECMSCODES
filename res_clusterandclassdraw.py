#!/usr/bin/env python
# -*- coding: utf-8 -*-
# @Author : Ji Hongchen
# @Email : jhca.xyt@163.com
# @Last version : 2021-12-27 23:37
# @Filename : pl_drawpilemap.py

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import cfg
import seaborn as sns

CFG = cfg.cfg()

ITEM = 'COMB'
#DRAW_COL = 'gender'
DRAW_COL = 'tissue_or_organ_of_origin'
#DRAW_COL = 'race'
#DRAW_COL = 'ajcc_pathologic_stage'
#DRAW_COL = 'ajcc_pathologic_t'
#DRAW_COL = 'ajcc_pathologic_n'
#DRAW_COL = 'ajcc_pathologic_m'
#DRAW_COL = 'tissue_or_organ_of_origin'


CLINICAL_FILE = CFG.datapath + 'clinical_TCGA_for_R.tsv'
CLASS_FILE = CFG.lassoresultpath + 'LASSO_CLASS_COMB.csv'

colors = CFG.colors


def DrawPileMap(dataframe, class_col, draw_col, draw_order='nan'):
    hue_class = list(set(dataframe[class_col]))
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
    plt.ylim(0, 1)
    plt.xticks(range(len(hue_class)), hue_class, rotation='vertical', size=20)
    plt.yticks([0, 0.25, 0.5, 0.75, 1],
               ['0%', '25%', '50%', '75%', '100%'],
               size=20)
    plt.legend(loc=[0, 1.1])
    return plt


if __name__ == '__main__':
    clinicaldf = pd.read_csv(CLINICAL_FILE, sep='\t', index_col=[0])
    print(clinicaldf.columns)
    classinfo = pd.read_csv(CLASS_FILE, index_col=[0])[[ITEM]]
    clinicaldf['define'] = 'a'
    for i in clinicaldf.index:
        try:
            clinicaldf['define'].loc[i] = classinfo[ITEM][i]
        except KeyError:
            pass
    clinicaldf.drop(clinicaldf[clinicaldf['define']
                    == 'a'].index, inplace=True)

    plt = DrawPileMap(clinicaldf, 'define', DRAW_COL)
    plt.tight_layout()
    plt.savefig(CFG.lassoresultpath + ITEM + '_' + DRAW_COL + '_pilemap.jpg')
    plt.close()


