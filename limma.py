#!/usr/bin/env python
# -*- coding: utf-8 -*-
# @Author : Ji Hongchen
# @Email : jhca.xyt@163.com
# @Last version : 2021-10-07 22:27
# @Filename : limmapreprocess.py

import pandas as pd
import numpy as np
import argparse
import os
import cfg
from DrawBoxPlots import _ExchangePatientID


CFG = cfg.cfg()
n_class = '2'

workspace = CFG.resultpath + 'SIG_RESULT/'
EXP_FILE = CFG.datapath + 'NORM_SELECTED_MRNA_SYMBOL.csv'
GSVA_FILE = CFG.datapath + 'gsva_result.csv'


def ChangeColNameOfGsvaResult(gsva_file):
    gf = pd.read_csv(gsva_file, index_col=[0])
    cols = gf.columns
    n_c = []
    for c in cols:
        if 'X' not in c:
            n_c.append(c)
        else:
            n_c.append(c.split('X')[1])
    l_c = []
    for c in n_c:
        l_c.append('-'.join(c.split('.')))
    print(l_c)
    gf.columns = l_c
    print('New frame:')
    print(gf)
    gf.to_csv(gsva_file)


def GenerateLimmaSeq(gsva_file, clas_file):
    gsva_result = pd.read_csv(gsva_file, index_col=[0])
    clas_result = pd.read_csv(clas_file, index_col=[0])
    class_frame = pd.DataFrame(columns=['class'])
    for i in gsva_result.columns:
        class_id = clas_result['finalclass'][i]
        class_frame.loc[i] = ['c_' + str(class_id)]
    return class_frame


if __name__ == '__main__':
    folderlist = os.listdir(workspace)
    n_cluster = [3, 4]
    for folder in folderlist:
        print('================================================')
        print('Folder: ', folder)
        gsva_result = workspace + folder + 'gsva_result.csv'
        for n in n_cluster:
            print('------------------------------------------------')
            print('n_cluster: ', n)
            classdict = np.load(
                workspace + folder + '/classdict_' + str(n) + '.npy', allow_pickle=True).item()
            cluster_result = pd.read_csv(
                workspace + folder + '/ClusterResult/cluster_' + str(n) + '.csv')
            for item in classdict:
                cluster_result['x'][cluster_result['x']
                                    == item] = classdict[item]
            print('>>> GSVA process...')
            os.system('Rscript limma.R ' +
                      GSVA_FILE + ' ' +
                      CFG.datapath + 'temp.csv ' +
                      str(n) + ' '
                      + workspace + folder + '/ClusterResult/gsva_limmaresult_' + str(n) + '.csv')
            for i in cluster_result.index:
                cluster_result['Unnamed: 0'][i] = _ExchangePatientID(
                    cluster_result['Unnamed: 0'][i])
            cluster_result.to_csv(CFG.datapath + 'temp.csv', index=False)
            print('>>> mRNA expression process...')
            os.system('Rscript limma.R ' +
                      EXP_FILE + ' ' +
                      CFG.datapath + 'temp.csv ' +
                      str(n) + ' '
                      + workspace + folder + '/ClusterResult/limmaresult_' + str(n) + '.csv')

    # ChangeColNameOfGsvaResult(GSVA_FILE)
    #rs = GenerateLimmaSeq(GSVA_FILE, CLAS_FILE)
    # rs.to_csv('limma_class.csv')
    #os.system('Rscript limma.R ' + GSVA_FILE + ' limma_class.csv ' + n_class)
    # os.remove('limma_class.csv')
