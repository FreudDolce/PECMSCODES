#!/usr/bin/env python
# -*- coding: utf-8 -*-
# @Author : Ji Hongchen
# @Email : jhca.xyt@163.com
# @Last version : 2021-10-28 10:33
# @Filename : SearchSurvivalGroup.py

import os
import pandas as pd
import numpy as np
import cfg

CFG = cfg.cfg()

SEARCHPATH = CFG.resultpath + 'RAND_RESULT/'


def GetClassForPatients(resultpath, n_cluster):
    resultfolderlist = os.listdir(resultpath)
    cols = ['case_id']
    for resultfolder in resultfolderlist:
        cols.append(resultfolder)
    resultframe = pd.DataFrame(columns=cols)
    for resultfolder in resultfolderlist:
        if os.path.isdir(resultpath + resultfolder) == True:
            resultfile = resultpath + resultfolder + \
                '/ClusterResult/With_class_cluster_' + str(n_cluster) + '.csv'
            df = pd.read_csv(resultfile, index_col=[0])
            df[resultfolder] = df['class']
            pdf = df[['case_id', resultfolder]]
            try:
                total_class_list = pd.merge(
                    total_class_list, pdf, how='inner', on='case_id')
            except UnboundLocalError:
                total_class_list = pdf
    total_class_list.to_csv(SEARCHPATH + str(n_cluster) +
                            '_merge_class.csv', index=False)
    return total_class_list


def _dicecoef(list_1, list_2):
    common_list = []
    for i in list_1:
        if i in list_2:
            common_list.append(i)
    dicecoef = len(common_list) * 2 / (len(list_1) + len(list_2))
    return dicecoef


def GetClassNuber(classfile):
    dataframe = pd.read_csv(classfile, index_col=[0])
    cols = list(dataframe.columns)
    for c in dataframe.columns:
        dataframe['temp'] = 0
        dataframe['temp'] = dataframe[c].str[-1]
        dataframe[c] = dataframe['temp']
    del dataframe['temp']
    dataframe['finalclass'] = 0
    dataframe[cols] = dataframe[cols].astype('int')
    for i in dataframe.index:
        dataframe['finalclass'].loc[i] = \
                dataframe[cols].loc[i].sum() / len(cols)
    dataframe['finalclass'][dataframe['finalclass'] > 1.8] = 2
    dataframe['finalclass'][dataframe['finalclass'] < 1.2] = 1
    dataframe.to_csv(SEARCHPATH + 'Prob_class.csv')
    dataframe.drop(dataframe[dataframe['finalclass'].isin([1.25, 1.5, 1.75])].index, inplace=True)
    dataframe.to_csv(SEARCHPATH + 'Cert_class.csv')


if __name__ == '__main__':
    #ss = GetClassForPatients(SEARCHPATH, 2)
    GetClassNuber(SEARCHPATH + '2_merge_class.csv')
