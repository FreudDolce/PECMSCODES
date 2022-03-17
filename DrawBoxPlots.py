#!/usr/bin/env python
# -*- coding: utf-8 -*-
# @Author : Ji Hongchen
# @Email : jhca.xyt@163.com
# @Last version : 2021-11-08 19:11
# @Filename : estimateanalysis.py

import pandas as pd
import numpy as np
import os
import matplotlib
import matplotlib.pyplot as plt
import seaborn as sns
import cfg
import argparse

CFG = cfg.cfg()

parser = argparse.ArgumentParser()
parser.add_argument('-n', help='Number of cluster')
args = parser.parse_args()

def _ExchangePatientID(patientID):
    if 'X' in patientID:
        newpid = patientID.split('X')[1]
        newpid = '-'.join(newpid.split('.'))
    else:
        newpid = '-'.join(patientID.split('.'))
    return newpid


def MergeClassFile(datafile, classfile, inplacedict=[], has_inplace=False, has_changeid=False):
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
    plt.figure(figsize=(0.4 * len(statitem) + 0.5, 7.5), dpi=300)
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
    ana_folder = CFG.resultpath + 'SIG_RESULT/'
    folders = os.listdir(ana_folder)
    savename = 'purity-'
    #items = ['B cells memory', 'Plasma cells', 'T cells CD8',
    #         'T cells CD4 memory resting', 'T cells CD4 memory activated',
    #         'T cells regulatory (Tregs)', 'NK cells activated',
    #         'Macrophages M1', 'Macrophages M2', 'Neutrophils']
    #items = ['StromalScore', 'ImmuneScore', 'ESTIMATEScore']
    items = ['TumorPurity']
    for folder in folders:
        nf = MergeClassFile(
            CFG.datapath + 'estimate_score.csv', #'estimate_score.csv',
            ana_folder + folder + '/ClusterResult/cluster_' + str(args.n) + '.csv',
            ana_folder + folder + '/classdict_' + str(args.n) + '.npy',
            has_inplace=True,
            has_changeid=False)
        pf = Platedata(nf)
        plt = DrawMutiVioPlot(dataframe=pf,
                              statcols=['hue', 'classname'],
                              valuecol='value',
                              statitem=items)
        plt.xticks(rotation=90)
        plt.tight_layout()
        plt.savefig(ana_folder + folder + '/' + savename + str(args.n) + '_vio.jpg')
        plt.close()
        plt = DrawMutiBoxPlot(dataframe=pf,
                              statcols=['hue', 'classname'],
                              valuecol='value',
                              statitem=items)
        plt.xticks(rotation=90)
        plt.tight_layout()
        plt.savefig(ana_folder + folder + '/' + savename + str(args.n) + '_box.jpg')
        plt.close()
