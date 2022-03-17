#!/usr/bin/env python
# -*- coding: utf-8 -*-
# @Author : Ji Hongchen
# @Email : jhca.xyt@163.com
# @Last version : 2021-12-08 16:06
# @Filename : difgeneanalysis.py

import pandas as pd
import numpy as np
import os
import cfg
from DrawBoxPlots import _ExchangePatientID
from scipy.stats import ttest_ind

CFG = cfg.cfg()


def GetGroupList(dataframe, classframe, log=True):
    if log == True:
        dataframe = np.log2(dataframe)
    classlist = list(set(classframe['x']))
    exp_dict = {}
    for c in classlist:
        sample_list = []
        sample_list_ori = list(
            classframe['Unnamed: 0'][classframe['x'] == c])
        for s in sample_list_ori:
            sample_list.append(_ExchangePatientID(s))
        exp_dict[c] = dataframe[sample_list]
    return exp_dict


def GeneExpTtest(df1, df2, compare_list, has_log=True):
    test_result = pd.DataFrame(
        columns=['avr_1', 'n_1', 'avr_2', 'n_2', 'log2fc', 'p', 'fdr'])
    for gene in compare_list:
        print('>>> gene: ', gene)
        genedata1 = df1.loc[gene]
        genedata2 = df2.loc[gene]
        avr_1 = genedata1.mean()
        avr_2 = genedata2.mean()
        n_1 = len(genedata1)
        n_2 = len(genedata2)
        log2fc = avr_2 - avr_1
        p = ttest_ind(genedata1, genedata2)[1]
        fdr = 'nan'
        test_result.loc[gene] = [
            avr_1, n_1, avr_2, n_2, log2fc, p, fdr]
    return test_result


if __name__ == '__main__':
    comb_dict = {'COL_H-COL-L': {'COL_H': ['COL', 'D'], 'COL_L': ['HA', 'L']},
                 'HA_H-HA-L': {'HA_H': ['HA', 'D'], 'HA_L': ['COL', 'L']},
                 'HA-COL': {'HA_ALONE': ['HA'], 'COL_ALONE': ['COL']}}
    exp_data = pd.read_csv(
        CFG.datapath + 'SELECTED_MRNA_SYMBOL.csv', index_col=[0])
    genelist = list(exp_data.index)
    workspace = CFG.resultpath + 'SIG_RESULT/'
    for folder in os.listdir(workspace):
        class_data = pd.read_csv(
            workspace + folder + '/ClusterResult/cluster_4.csv')
        name_dict = np.load(
            workspace + folder + '/classdict_4.npy',
            allow_pickle=True).item()
        for n in name_dict:
            class_data['x'][class_data['x'] == n] = name_dict[n]
        for test_item in comb_dict:
            combclassdata = class_data.copy()
            print('...', test_item)
            for item in comb_dict[test_item]:
                combclassdata['x'][combclassdata['x'].isin(
                    comb_dict[test_item][item])] = item
            combclassdata = combclassdata[combclassdata['x'].isin(
                comb_dict[test_item])]
            exp_dict = GetGroupList(exp_data, combclassdata, log=True)
            exp_test_list = GeneExpTtest(exp_dict[list(exp_dict)[0]],
                                         exp_dict[list(exp_dict)[1]],
                                         compare_list=genelist,
                                         has_log=True)
            exp_test_list.to_csv(workspace + folder + '/ClusterResult/' + \
                    test_item + '_exp_dif.csv')
