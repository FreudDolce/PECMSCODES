#!/usr/bin/env python
# -*- coding: utf-8 -*-
# @Author : Ji Hongchen
# @Email : jhca.xyt@163.com
# @Last version : 2022-01-11 15:37
# @Filename : res_combine_cpt_tcga_data.py

import pandas as pd
import numpy as np

TYPE = 'TYPE1'
coef = '2'
inplace = True

TCGA_FILE = pd.read_csv(
    '~/Documents/MANU/BI/PANC_ECM/Manuscript/Data/Clinica_TCGA_' +
    TYPE + '_COEF_' + coef + '.csv',
    index_col=[0])
CPT_FILE = pd.read_csv(
    '~/Documents/MANU/BI/PANC_ECM/Manuscript/Data/Clinica_CPT_' +
    TYPE + '_COEF_' + coef + '.csv',
    index_col=[0])

COMB_FILE = pd.concat([TCGA_FILE, CPT_FILE], join='inner')
COMB_FILE.to_csv('~/Documents/MANU/BI/PANC_ECM/Manuscript/Data/Clinica_COMB_' +
                 TYPE + '_COEF_' + coef + '.csv')
