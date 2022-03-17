#!/usr/bin/env python
# -*- coding: utf-8 -*-
# @Author : Ji Hongchen
# @Email : jhca.xyt@163.com
# @Last version : 2021-10-01 23:42
# @Filename : searchpathway.py

import pandas as pd
import numpy as np
import os
import argparse
import cfg

CFG = cfg.cfg()

parser = argparse.ArgumentParser()
parser.add_argument('-p', type=str, help='Pathway you want to serch')
parser.add_argument('-e', type=str, default='CCXXCC_XXCCXX', help='Pathway you want to exculded.')
parser.add_argument('-s', type=str, help='if 0, donot save, else, s is the saving folder name')
args = parser.parse_args()

ENZ2NAME = np.load(CFG.GENE_ID_TRANS_PATH + 'enz2name.npy',
                   allow_pickle=True).item()
PATH_PATH = CFG.SELECTED_PATHWAY_PATH


def GeneratePathDict(gmt_file):
    gmt_info = pd.read_csv(gmt_file, names=['p'])
    gmt_info['pathwayname'] = 'a'
    gmt_info['pathwayname'] = gmt_info['p'].str.split('\t', expand=True)[0]
    pathwaylist = list(gmt_info['pathwayname'])
    np.save('pathwaylist.npy', pathwaylist)


def SearchPathway(searchname, excludename):
    nameplit = []
    plist = np.load('/Users/freud/Documents/MANU/BI/PANC_ECM/DATA/gmt/pathwaylist.npy')
    for p in plist:
        if searchname in p:
            if excludename not in p:
                nameplit.append(p)
    return nameplit


def CalcMoleDict(pathwayclass):
    genelist = []
    pathway_list = pd.read_csv(
        PATH_PATH + pathwayclass + '_search_result.csv')['pathway']
    for i in pathway_list:
        for j in pathway_filelist:
            try:
                genelist.append(j[i])
            except KeyError:
                pass
    return genelist


if __name__ == '__main__':
    #gl = CalcMoleDict('ATP')
    # GeneratePathDict(CFG.ORI_GMT)
    # """
    searched_pathway = SearchPathway(args.p, args.e)
    searched_frame = pd.DataFrame(searched_pathway, columns=['pathway'])
    print(searched_frame)
    if args.s != '0':
        if os.path.exists(CFG.SEARCHED_PATHWAY + args.s) == False:
            os.mkdir(CFG.SEARCHED_PATHWAY + args.s)
        searched_frame.to_csv(CFG.SEARCHED_PATHWAY + args.s + '/' + args.p + '_search_result.csv',
                              index=False)
    # """
