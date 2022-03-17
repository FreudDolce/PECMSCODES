#!/usr/bin/env python
# -*- coding: utf-8 -*-
# @Author : Ji Hongchen
# @Email : jhca.xyt@163.com
# @Last version : 2021-11-11 23:46
# @Filename : FALL.py

import pandas as pd
import numpy as np
import os
from DrawBoxPlots import _ExchangePatientID
import cfg

CFG = cfg.cfg()

workspace = CFG.resultpath + 'SIG_RESULT/'

clinicalannotation = CFG.datapath + 'clinical_with_history.tsv'


def MergeclassResult(clinicalfile, classinfo, classdict='none', changeclass=False):
    if changeclass == True:
        classname = np.load(classdict, allow_pickle=True).item()
    clinicalinfo = pd.read_csv(clinicalfile, sep='\t')
    classinfo['case_id'] = 'a'
    classinfo['define'] = 'a'
    for i in classinfo.index:
        classinfo['case_id'][i] = _ExchangePatientID(
            classinfo['Unnamed: 0'][i])
    if changeclass == True:
        for j in classname:
            classinfo['define'][classinfo['x'] == j] = classname[j]
    else:
        pass
    withclassinfo = pd.merge(
        clinicalinfo, classinfo[['case_id', 'define']], on='case_id')
    return withclassinfo


if __name__ == '__main__':
    """
    classinfo = pd.read_csv(CFG.lassoresultpath +
                            'LASSO_CLASS_COMB.csv', index_col=[0])
    classinfo = classinfo.reset_index()
    item = 'COL'
    classinfo.rename(columns={'index': 'case_id'}, inplace=True)
    cinfo = classinfo[['Unnamed: 0', item]]
    withclassframe = MergeclassResult(
        clinicalfile=clinicalannotation,
        classinfo=cinfo,
    )
    del withclassframe['define']
    withclassframe.rename(columns={item: 'define'}, inplace=True)
    print(withclassframe)
    """
    n_cluster = [3, 4]
    folderlist = os.listdir(workspace)
    for folder in folderlist:
        print('============================================')
        print('Using result: ', folder)
        for n in n_cluster:
            print('-------------------------------------------')
            print('number of class: ', n)
            cfile = workspace + folder + \
                '/ClusterResult/cluster_' + str(n) + '.csv'
            cinfo = pd.read_csv(cfile)
            cdict = workspace + folder + '/classdict_' + str(n) + '.npy'
            withclassframe = MergeclassResult(
                clinicalfile=clinicalannotation,
                classinfo=cinfo,
                classdict=cdict,
                changeclass=True)
            print(withclassframe)
            withclassframe.to_csv(
                CFG.datapath + 'temp.tsv', sep='\t', index=False)
            maffile = 'SomaticMutation/TCGA.PAAD.mutect.fea333b5-78e0-43c8-bf76-4c78dd3fac92.DR-10.0.somatic.maf.gz'
            figsavepath = workspace + folder + '/ClusterResult/mutect_' + \
                str(n) + '.jpg'
            os.system('Rscript ' + CFG.codepath +
                      'FALL.R ' + maffile + ' temp.tsv ' + figsavepath)
            os.remove(CFG.datapath + 'temp.tsv')
