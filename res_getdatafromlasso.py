#!/usr/bin/env python
# -*- coding: utf-8 -*-
# @Author : Ji Hongchen
# @Email : jhca.xyt@163.com
# @Last version : 2022-01-04 01:16
# @Filename : res_getdatafromlasso.py

import pandas as pd
import numpy as np

clicinfo = pd.read_csv('~/Documents/MANU/BI/PANC_ECM/Manuscript/Data/TCGA_ClinicalwithCluster.csv',
                       index_col=[0])
cols = list(clicinfo)
cols[-1] = 'ori_comb'
clicinfo.columns = cols

exp_info = pd.read_csv('~/Documents/MANU/BI/PANC_ECM/DATA/SELECTED_MRNA_SYMBOL.csv',
                       index_col=[0])
exp_info = pd.DataFrame(np.array(exp_info).T,
                        index=exp_info.columns, columns=exp_info.index)

genelist = list(pd.read_csv('~/Documents/MANU/BI/PANC_ECM/Manuscript/Data/limma_COMB_sig.csv',
                            index_col=[0]).index)
exp_info = exp_info[genelist]
exp_info = np.log2(exp_info)

exp_info = exp_info.join(clicinfo[['ori_comb']])
exp_info_2 = exp_info[exp_info['ori_comb'].isin(['L', 'D'])]
exp_info_2['ori_comb'][exp_info_2['ori_comb'] == 'D'] = 1
exp_info_2['ori_comb'][exp_info_2['ori_comb'] == 'L'] = 0
print(exp_info_2)
exp_info_2.to_csv('~/Documents/MANU/BI/PANC_ECM/Manuscript/Data/comb_for_limma_norm.csv')
