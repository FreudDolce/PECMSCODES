#!/usr/bin/env python
# -*- coding: utf-8 -*-
# @Author : Ji Hongchen
# @Email : jhca.xyt@163.com
# @Last version : 2021-12-25 12:35
# @Filename : pl_fall.py

import pandas as pd
import numpy as np
import cfg
import os

CFG = cfg.cfg()
ITEM = 'COL'

maffile = 'SomaticMutation/cpt_comb.maf'
clinicalannotation = CFG.datapath + 'CPT_clinical.csv'
clinicalinfo = pd.read_csv(clinicalannotation)
classinfo = pd.read_csv(CFG.lassoresultpath + 'LASSO_CPT_COMB.csv')
classinfo.rename(columns={'Unnamed: 0': 'case_id'}, inplace=True)
classinfo.rename(columns={ITEM: 'define'}, inplace=True)
df = pd.merge(clinicalinfo, classinfo[['case_id', 'define']], on='case_id')
df.to_csv(CFG.datapath + 'temp.tsv', sep='\t', index=False)

figsavepath = CFG.lassoresultpath + ITEM + '_cpt_mutect_fall.jpg'
os.system('Rscript ' + CFG.codepath +
          'FALL.R ' + maffile + ' temp.tsv ' + figsavepath)
#os.remove(CFG.datapath + 'temp.tsv')
