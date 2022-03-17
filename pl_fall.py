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
ITEM = 'COMB'

clinicalannotation = CFG.datapath + 'clinical_with_history.tsv'
clinicalinfo = pd.read_csv(clinicalannotation, sep='\t')
classinfo = pd.read_csv(CFG.lassoresultpath + 'LASSO_CLASS_COMB.csv')
classinfo.rename(columns={'Unnamed: 0': 'case_id'}, inplace=True)
classinfo.rename(columns={ITEM: 'define'}, inplace=True)
df = pd.merge(clinicalinfo, classinfo[['case_id', 'define']], on='case_id')
df.to_csv(CFG.datapath + 'temp.tsv', sep='\t', index=False)

figsavepath = CFG.lassoresultpath + ITEM + '_mutect_fall.jpg'
maffile = 'SomaticMutation/TCGA.PAAD.mutect.fea333b5-78e0-43c8-bf76-4c78dd3fac92.DR-10.0.somatic.maf.gz'
os.system('Rscript ' + CFG.codepath +
          'FALL.R ' + maffile + ' temp.tsv ' + figsavepath)
os.remove(CFG.datapath + 'temp.tsv')
