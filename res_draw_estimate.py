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
    plt.figure(figsize=(1.5 * len(statitem) + 0.5, 5.5), dpi=300)
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
    #plt.figure(figsize=(0.4 * len(statitem) + 0.5, 7.5), dpi=300)
    sns.boxplot(x=statcols[0],
                y=valuecol,
                hue=statcols[1],
                data=drawframe,
                linewidth=0.5,
                width=0.7,
                fliersize=1.0
                )
    plt.legend(loc=8, bbox_to_anchor=[0, 1])
    return plt


if __name__ == '__main__':
    classfile = pd.read_csv(
        '~/Documents/MANU/BI/PANC_ECM/Manuscript/Data/' + tyep + '_post_lasso_coef_' + coef + '.csv')
    classfile = classfile[['Unnamed: 0', 'lasso_class']]
    nf = MergeClassFile(
        CFG.datapath + 'estimate_score_for_plot.csv',
        classframe=classfile,
        inplacedict=[],
        has_inplace=False,
        has_changeid=False)
    #items = ['TumorPurity']
    items = ['StromalScore', 'ImmuneScore', 'ESTIMATEScore']
    pf = Platedata(nf)
    figure = plt.figure(figsize=(3, 8), dpi=300)
    plt = DrawMutiBoxPlot(dataframe=pf,
                          statcols=['hue', 'classname'],
                          valuecol='value',
                          statitem=items)
    plt.xticks(rotation=90)
    plt.tight_layout()
    plt.savefig('/Users/freud/Documents/MANU/BI/PANC_ECM/Manuscript/Figure/SubFigure/' +
                tyep + '_coef_' + coef + '_' + items[0] + '.tiff')
    plt.close()
