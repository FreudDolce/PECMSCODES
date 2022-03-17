#!/usr/bin/env python
# -*- coding: utf-8 -*-
# @Author : Ji Hongchen
# @Email : jhca.xyt@163.com
# @Last version : 2021-12-06 00:13
# @Filename : bi_genelimma.py

import pandas as pd
import numpy as np
import os
import cfg
from DrawBoxPlots import _ExchangePatientID

CFG = cfg.cfg()

if __name__ == '__main__':
    exp_data = CFG.datapath + 'SELECTED_MRNA_SYMBOL.csv'
    workspace = CFG.resultpath + 'SIG_RESULT/'
    folderlist = os.listdir(workspace)
    n_cluster = [4]
    for folder in folderlist:
        print('Using folder: >>> ', folder)
        for n in n_cluster:
            print('Number of classes: ', n)
            class_dict = np.load(
                workspace + folder + '/classdict_' + str(n) + '.npy',
                allow_pickle=True).item()
            class_file = workspace + folder + \
                '/ClusterResult/' + 'cluster_' + str(n) + '.csv'
            class_info = pd.read_csv(class_file)
            for ind in class_info.index:
                class_info['Unnamed: 0'][ind] = \
                    _ExchangePatientID(class_info['Unnamed: 0'][ind])
            for i in class_dict:
                class_info['x'][class_info['x'] ==i] = class_dict[i]
            class_info = class_info[class_info['x'].isin(['HA', 'COL'])]
            print(class_info.shape)
            # COL_vs_HA
            mrna_data = pd.read_csv(exp_data)
            cols = ['ENS']
            cols.extend(list(class_info['Unnamed: 0']))
            mrna_data = mrna_data[cols]
            print(mrna_data.shape)
            mrna_data.to_csv('mrna_temp.csv', index=False)
#            class_info['x'][class_info['x'].isin(['HA', 'D'])] = 'HA_H'
#            class_info['x'][class_info['x'].isin(['COL', 'L'])] = 'HA_L'
            print(class_info)
            class_info.to_csv('classlist.csv', index=False)
            filename = workspace + folder + '/ClusterResult/'\
                    'COL_vs_HA_limma.csv'
            os.system('Rscript ' + CFG.codepath + 'limma.R mrna_temp.csv classlist.csv ' +\
                    '2' + ' ' + filename)
            os.remove('classlist.csv')
            os.remove('mrna_temp.csv')
