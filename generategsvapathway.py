#!/usr/bin/env python
# -*- coding: utf-8 -*-
# @Author : Ji Hongchen
# @Email : jhca.xyt@163.com
# @Last version : 2021-10-02 19:35
# @Filename : generategsvapathway.py

import pandas as pd
import numpy as np
import os
import argparse
import random
import cfg

CFG = cfg.cfg()

SELECT_PATH = CFG.datapath + 'SelectedGmt/'
KEGFILE = CFG.worksapce + 'DATA/gmt/msigdb_v7.4/msigdb_v7.4_GMTs/msigdb.v7.4.entrez.gmt'

parser = argparse.ArgumentParser()
parser.add_argument('-s', help='Save name of gmk file')
args = parser.parse_args()

filelist = os.listdir(SELECT_PATH)
kwargs = {'sep': '\t'}


def JointPathway():
    pathway_list = []
    for f in filelist:
        pathways = np.array(pd.read_csv(SELECT_PATH + f)['pathway'])
        pathway_list.extend(list(pathways))

    selectedframe = pd.DataFrame(columns=['sp'])
    for i in range(len(keginfo)):
        for j in pathway_list:
            if j in keginfo['p'][i]:
                print(keginfo['p'][i])
                selectedframe.append(keginfo.loc[i])


def Getpathway(pathwaypath):
    pathway_list = []
    fl = os.listdir(pathwaypath)
    for f in fl:
        if 'limma' in f:
            ss = pd.read_csv(pathwaypath + f, index_col=[0])
            pathways = list(ss.index)
            for p in pathways:
                pathway_list.append(p.split(',')[1])
        elif 'search_result' in f:
            ss = pd.read_csv(pathwaypath + f)
            pathways = list(ss['pathway'])
            for p in pathways:
                pathway_list.append(p)
        else:
            print('Error file: ', f)
    print('Not merged: ', len(pathway_list), ' pathways..')
    pathway_list = list(set(pathway_list))
    return pathway_list


def GetGmtFile(pathwaylist, savename):
    print(KEGFILE)
    pathinfo = pd.read_csv(KEGFILE, names=['p'])
    pathinfo['tem'] = 'a'
    pathinfo['tem'] = pathinfo['p'].str.split('\t', expand=True)[0]
    newpathwayframe = pathinfo[pathinfo['tem'].isin(pathwaylist)]
    newpathwayframe = newpathwayframe[['p']][0:]
    print(newpathwayframe)
    print('Selected ', len(newpathwayframe), ' pathways..')
    newpathwayframe.to_csv(SELECT_PATH + savename + '.gmt', index=False, header=0)


def MergeSearchedPathway(folder_path):
    for folder in os.listdir(folder_path):
        for f in os.listdir(folder_path + '/' + folder):
            finfo = pd.read_csv(folder_path + '/' + folder + '/' + f)
            try:
                folderinfo = folderinfo.append(finfo)
            except NameError:
                folderinfo = finfo
        folderinfo.reset_index(inplace=True, drop=True)
        folderinfo.to_csv(CFG.MERGED_PATH + folder + '.csv', index=False)
        del folderinfo


def RandomSelectPathway(path_path):
    flist = os.listdir(path_path)
    number_of_cluster = random.randint(1, 3)
    pathway_per_cluster = int(CFG.PATH_SELECT_PER_ITER / number_of_cluster)
    selected_file = random.sample(flist, number_of_cluster)
    print('Selected cluster: ', selected_file, ', ', number_of_cluster, ' in total.')
    for f in selected_file:
        finfo = pd.read_csv(path_path + '/' + f)
        selected_finfo = finfo.loc[random.sample(range(len(finfo)), pathway_per_cluster)]
        try:
            all_selected_finfo = all_selected_finfo.append(selected_finfo)
        except NameError:
            all_selected_finfo = selected_finfo
    all_selected_finfo.reset_index(inplace=True, drop=True)
    all_selected_finfo_list = list(set(all_selected_finfo['pathway']))
    return all_selected_finfo_list



if __name__ == '__main__':
    # MergeSearchedPathway(CFG.SEARCHED_PATHWAY)
    #selected_pathway = RandomSelectPathway(CFG.MERGED_PATH)
    gmts_path = CFG.datapath + 'Searched_pathway/PD/'
    fl = os.listdir(gmts_path)
    for f in fl:
        finfo = pd.read_csv(gmts_path + f)
        try:
            tinfo = tinfo.append(finfo)
        except NameError:
            tinfo = finfo
    plist = list(set(tinfo['pathway']))
    # pl = Getpathway(SELECT_PATH)
    GetGmtFile(plist, args.s)
