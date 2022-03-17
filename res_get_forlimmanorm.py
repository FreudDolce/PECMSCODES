#!/usr/bin/env python
# -*- coding: utf-8 -*-
# @Author : Ji Hongchen
# @Email : jhca.xyt@163.com
# @Last version : 2022-01-08 23:00
# @Filename : res_get_forlimmanorm.py

import pandas as pd
import numpy as np

item = 'TYPE2'

genelist = pd.read_csv('~/Documents/MANU/BI/PANC_ECM/Manuscript/Data/' +
                       item + '_limma_sig.csv',
                       index_col=[0])
#genelist = pd.read_csv('~/Documents/MANU/BI/PANC_ECM/Manuscript/Data/SingleCoxSig.csv',
#                       index_col=[0])

clic = pd.read_csv(
    '~/Documents/MANU/BI/PANC_ECM/Manuscript/Data/TCGA_ClinicalwithCluster.csv',
    index_col=[0])[['ori_class']]

genelist = list(genelist.index)
exp = pd.read_csv(
    '~/Documents/MANU/BI/PANC_ECM/DATA/SELECTED_MRNA_SYMBOL.csv', index_col=[0])

r_exp = pd.DataFrame(np.array(exp).T, index=exp.columns, columns=exp.index)
r_exp = r_exp[genelist]
r_exp = np.log2(r_exp)
r_exp = r_exp.join(clic)
r_exp['ori_class'][r_exp['ori_class'] == 'D'] = 1
r_exp['ori_class'][r_exp['ori_class'] == 'HA'] = 1
r_exp['ori_class'][r_exp['ori_class'] != 1] = 0
print(r_exp)

r_exp.to_csv('~/Documents/MANU/BI/PANC_ECM/Manuscript/Data/' +
             item + '_for_limma_norm.csv')
