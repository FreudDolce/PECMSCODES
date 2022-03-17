#!/usr/bin/env python
# -*- coding: utf-8 -*-
# @Author : Ji Hongchen
# @Email : jhca.xyt@163.com
# @Last version : 2021-12-17 22:36
# @Filename : checkdifgene.py

import pandas as pd
import numpy as np
import cfg
import os

CFG = cfg.cfg()
P_CUT = 0.00001
FC_CUT = 1.3

workspace = CFG.resultpath + \
    'SIG_RESULT/0.0021-str-comb_20211125/ClusterResult/'
DIF_FILE = workspace + 'HA_H-HA-L_exp_dif.csv'

def getselectedgenes(file_name, p_cut, fc_cut):
    dif_info = pd.read_csv(file_name, index_col=[0])
    dif_info['abs'] = dif_info['log2fc'].abs()
    dif_info.sort_values(by='abs', ascending=False, inplace=True)
    dif_info = dif_info[dif_info['abs'] > fc_cut]
    dif_info = dif_info[dif_info['fdr'] < p_cut]
    return dif_info

if __name__ == '__main__':
    di = getselectedgenes(DIF_FILE, P_CUT, FC_CUT)
    print(di)
    print(di.shape)
    di.to_csv(workspace + 'HA_sig_dif.csv')
