#!/usr/bin/env python
# -*- coding: utf-8 -*-
# @Author : Ji Hongchen
# @Email : jhca.xyt@163.com
# @Last version : 2021-10-10 22:10
# @Filename : RandomSearchpathway.py

import os
import argparse
import pandas as pd
import cfg
import shutil
import time
import random
import numpy as np

CFG = cfg.cfg()


SELECT_PATH = CFG.datapath + 'SelectedGmt/'
KEGFILE = CFG.datapath + 'gmt/msigdb_v7.4/msigdb_v7.4_GMTs/msigdb.v7.4.entrez.gmt'


def RandomSelectPathway(path_path):
    flist = os.listdir(path_path)
    number_of_cluster = random.randint(1, 3)
    pathway_per_cluster = int(CFG.PATH_SELECT_PER_ITER / number_of_cluster)
    selected_file = random.sample(flist, number_of_cluster)
    print('Selected cluster: ', selected_file,
          ', ', number_of_cluster, ' in total.')
    for f in selected_file:
        finfo = pd.read_csv(path_path + '/' + f)
        selected_finfo = finfo.loc[random.sample(
            range(len(finfo)), pathway_per_cluster)]
        try:
            all_selected_finfo = all_selected_finfo.append(selected_finfo)
        except NameError:
            all_selected_finfo = selected_finfo
    all_selected_finfo.reset_index(inplace=True, drop=True)
    all_selected_finfo_list = list(set(all_selected_finfo['pathway']))
    pathway_info = []
    for s in selected_file:
        pathway_info.append(s.split('.')[0])
    return ('_'.join(pathway_info), all_selected_finfo_list)


def GetGmtFile(pathwaylist, savename):
    pathinfo = pd.read_csv(KEGFILE, names=['p'])
    pathinfo['tem'] = 'a'
    pathinfo['tem'] = pathinfo['p'].str.split('\t', expand=True)[0]
    newpathwayframe = pathinfo[pathinfo['tem'].isin(pathwaylist)]
    newpathwayframe = newpathwayframe[['p']][0:]
    print('Selected ', len(newpathwayframe), ' pathways..')
    newpathwayframe.to_csv(SELECT_PATH + savename +
                           '.gmt', index=False, header=0)


def GenerateRandomResult(gmtfile, selected_cluster, time):
    # Make a dir: ClusterResult:
    if os.path.exists(CFG.CLUSTER_RESULT_PATH) == True:
        shutil.rmtree(CFG.CLUSTER_RESULT_PATH)
    os.mkdir(CFG.CLUSTER_RESULT_PATH)

    # Generate gsva sorce .csv data, cluster result and gsva figure, in CFG.resultpath
    os.system(CFG.codepath + 'gsva.py -g ' + gmtfile + '.gmt')

    # Merger clinical data into gsva cluster result, in CFG.CLUSTER_RESULT_PATH
    os.system(CFG.codepath + 'SurvivalAnalysis.py')

    # Calculate suvival P value Chi2 value, store in CFG.CLUSTER_RESULT_PATH /Survival_p.csv
    for i in range(CFG.MIN_CLASS, CFG.MAX_CLASS + 1):
        os.system('Rscript ' + CFG.codepath + 'SurvivalAnalysis.R ' +
                  CFG.CLUSTER_RESULT_PATH + 'With_class_cluster_' + str(i) + '.csv')
        survival_stat = pd.read_csv(CFG.codepath + 'survival_p.csv')
        try:
            survival_frame = survival_frame.append(
                survival_stat.loc[0], ignore_index=True)
        except NameError:
            survival_frame = survival_stat
        os.remove(CFG.codepath + 'survival_p.csv')
    survival_frame.drop('Unnamed: 0', axis=1)
    min_survival_p = round(survival_frame['V2'].min(), 4)
    print('Mini p value of survival analysis:', min_survival_p)
    survival_frame.to_csv(CFG.resultpath + 'Survival_p.csv',
                          index=range(CFG.MIN_CLASS + 1, CFG.MAX_CLASS + 1))

    os.system(CFG.codepath + 'CIBERSORT.py')

    CLEAR_PATH = CFG.resultpath
    if min_survival_p < 0.05:
        ORDER_PATH = CFG.resultpath + 'SIG_RESULT/'
    else:
        ORDER_PATH = CFG.resultpath + 'NSIG_RESULT/'

    fl = os.listdir(CLEAR_PATH)

    RESULT_PATH = ORDER_PATH + \
        str(min_survival_p) + '-' + selected_cluster + '_' + time + '/'
    os.mkdir(RESULT_PATH)
    for f in fl:
        if os.path.isfile(CLEAR_PATH + f) == True:
            shutil.move(CLEAR_PATH + f, RESULT_PATH + f)

    shutil.move(CLEAR_PATH + 'ClusterResult',
                RESULT_PATH + 'ClusterResult')


# Main programmed:
if __name__ == '__main__':
    while True:
        # Generate random selected gmt pathways
        t = time.localtime()
        now = str(t.tm_year) + '_' + str(t.tm_mon) + '_' + \
            str(t.tm_mday) + '_' + str(t.tm_hour) + '_' + str(t.tm_min)
        print('Now time: ', now)
        selected_file, selected_pathway = RandomSelectPathway(CFG.MERGED_PATH)
        GetGmtFile(selected_pathway, now)

        # Generate result, with folder name "now":
        GenerateRandomResult(now, selected_file, now)
