#!/usr/bin/env python
# -*- coding: utf-8 -*-
# @Author : Ji Hongchen
# @Email : jhca.xyt@163.com
# @Last version : 2021-11-28 18:57
# @Filename : genelimma.py

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
    n_cluster = [3, 4]
    for folder in folderlist:
        print('Using folder: >>> ', folder)
        for n in n_cluster:
            print('Number of classes: ', n)
            class_dict = np.load(
                workspace + folder + '/classdict_' + str(n) + '.npy',
                allow_pickle=True).item()
            class_file = workspace + folder + \
                '/ClusterResult/' + 'cluster_' + str(n) + '.csv'
            for i in class_dict:
                class_info = pd.read_csv(class_file)
                class_info['x'][class_info['x'] != i] = 'ELSE'
                class_info['x'][class_info['x'] != 'ELSE'] = class_dict[i]
                for ind in class_info.index:
                    class_info['Unnamed: 0'][ind] = \
                        _ExchangePatientID(class_info['Unnamed: 0'][ind])
                class_info.to_csv('classlist.csv', index=False)
                filename = workspace + folder + '/ClusterResult/'\
                        'Cluster_' + str(n) + '_' + class_dict[i] + '_vs_ELSE.csv'
                os.system('Rscript ' + CFG.codepath + 'limma.R ' + \
                        CFG.datapath + 'SELECTED_MRNA_SYMBOL.csv classlist.csv ' +\
                        '2' + ' ' + filename)
