#!/usr/bin/env python
# -*- coding: utf-8 -*-
# @Author : Ji Hongchen
# @Email : jhca.xyt@163.com
# @Last version : 2022-01-10 20:53
# @Filename : res_draw_cibersort.py

import pandas as pd
import numpy as np
import os
import matplotlib
import matplotlib.pyplot as plt
import seaborn as sns
import cfg
import argparse

CFG = cfg.cfg()
tyep = 'TYPE2'
coef = '2'


def _ExchangePatientID(patientID):
    if 'X' in patientID:
        newpid = patientID.split('X')[1]
        newpid = '-'.join(newpid.split('.'))
    else:
        newpid = '-'.join(patientID.split('.'))
    return newpid


def MergeClassFile(datafile, classframe, inplacedict=[], has_inplace=False, has_changeid=False):
    dataframe = pd.read_csv(datafile)
    classframe.columns = ['case_id', 'classname']
    if has_changeid == True:
        if '.' in classframe['case_id'][2]:
            for i in classframe.index:
                classframe['case_id'][i] = _ExchangePatientID(
                    classframe['case_id'][i])
    if has_inplace == True:
        classdict = np.load(inplacedict, allow_pickle=True).item()
        for cn in classdict:
            classframe['classname'][classframe['classname']
                                    == cn] = classdict[cn]
    neoframe = pd.merge(dataframe, classframe, on='case_id', how='inner')
    return neoframe


def Platedata(dataframe, store='classname'):
    cols = list(dataframe.columns)
    cols.remove(store)
    cols.remove('case_id')
    plateframe = pd.DataFrame(columns=['case_id', 'value', 'hue', 'classname'])
    for i in dataframe.index:
        for c in cols:
            plateframe = plateframe.append({'case_id': dataframe['case_id'][i],
                                            'value': dataframe[c][i],
                                            'hue': c,
                                            'classname': dataframe['classname'][i]},
                                           ignore_index=True)
    return plateframe


def DrawMutiVioPlot(dataframe, statcols, valuecol, statitem):
    drawframe = dataframe[dataframe['hue'].isin(statitem)]
    print(drawframe)
    sns.violinplot(x=statcols[0],
                   y=valuecol,
                   hue=statcols[1],
                   data=drawframe,
                   inner='point',
                   linewidth=0.5,
                   width=0.9
                   )
    plt.legend(loc=8, bbox_to_anchor=[0, 1])
    return plt


def DrawMutiBoxPlot(dataframe, statcols, valuecol, statitem):
    drawframe = dataframe[dataframe['hue'].isin(statitem)]
    sns.boxplot(x=statcols[0],
                y=valuecol,
                hue=statcols[1],
                data=drawframe,
                linewidth=0.5,
                width=0.8,
                fliersize=1.0
                )
    plt.legend(loc=8, bbox_to_anchor=[0, 1])
    return plt


if __name__ == '__main__':
    classfile = pd.read_csv(
        '~/Documents/MANU/BI/PANC_ECM/Manuscript/Data/' + tyep + '_post_lasso_coef_' + coef + '.csv')
    classfile = classfile[['Unnamed: 0', 'lasso_class']]
    nf = MergeClassFile(
        CFG.datapath + 'CIBERSORT_RESULT_for_plot.csv',
        classframe=classfile,
        inplacedict=[],
        has_inplace=False,
        has_changeid=False)
    print(nf.columns)
    pf = Platedata(nf)
    items = ['B cells naive', 'B cells memory', 'Plasma cells',
             'T cells CD8', 'T cells CD4 naive', 'T cells CD4 memory resting',
             'T cells CD4 memory activated', 'T cells follicular helper',
             'T cells regulatory (Tregs)', 'T cells gamma delta', 'NK cells resting',
             'NK cells activated', 'Monocytes', 'Macrophages M0', 'Macrophages M1',
             'Macrophages M2', 'Dendritic cells resting',
             'Dendritic cells activated', 'Mast cells resting',
             'Mast cells activated', 'Eosinophils', 'Neutrophils']
    plt.figure(figsize=(8, 8), dpi=300)
    plt = DrawMutiBoxPlot(dataframe=pf,
                          statcols=['hue', 'classname'],
                          valuecol='value',
                          statitem=items)
    bwith = 1.5  # 边框宽度设置为2
    ax = plt.gca()  # 获取边框
    ax.spines['bottom'].set_linewidth(bwith)
    ax.spines['left'].set_linewidth(bwith)
    ax.spines['top'].set_linewidth(0)
    ax.spines['right'].set_linewidth(0)
    plt.tick_params(width=1.5, length=5.5)
    plt.xlabel(' ')
    plt.ylabel(' ')
    plt.xticks(fontproperties='Arial', rotation=90, fontsize=20)
    plt.tight_layout()
    plt.yticks(fontproperties='Arial', size=20)
    plt.gcf().subplots_adjust(left=0.1, right=0.9, bottom=0.7)
    # plt.tight_layout()
    plt.savefig(
        '/Users/freud/Documents/MANU/BI/PANC_ECM/Manuscript/Figure/SubFigure/TYPE2_cibersort_box.tiff',
        transparent=True)
