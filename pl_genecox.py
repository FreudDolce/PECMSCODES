#!/usr/bin/env python
# -*- coding: utf-8 -*-
# @Author : Ji Hongchen
# @Email : jhca.xyt@163.com
# @Last version : 2021-12-28 17:21
# @Filename : pl_genecox.py

import pandas as pd
import numpy as np
import cfg
import matplotlib.pyplot as plt
from lifelines import CoxPHFitter

CFG = cfg.cfg()

SURV_FILE = '~/Documents/MANU/BI/PANC_ECM/DATA/TCGA_for_cox.csv'
EXP_FILE = CFG.datapath + 'SELECTED_MRNA_SYMBOL.csv'
expinfo = pd.read_csv(EXP_FILE, index_col=[0])
geneinfo = pd.DataFrame(np.array(expinfo).T,
                        index=expinfo.columns,
                        columns=expinfo.index)
survinfo = pd.read_csv(SURV_FILE)
survinfo = survinfo[['case_id', 'days_to_death', 'vital_status']]
survinfo.set_index('case_id', drop=True, inplace=True)

#col_genes = pd.read_csv('~/Documents/MANU/BI/PANC_ECM/Manuscript/Data/TYPE2_lasso_coef1.csv',
#                        index_col=[0])
#genelist = list(col_genes.index)
ha_genes = pd.read_csv('~/Documents/MANU/BI/PANC_ECM/Manuscript/Data/TYPE2_lasso_coef2.csv',
                       index_col=[0])
genelist = list(ha_genes.index)
#genelist.extend(list(ha_genes.index))
genelist = list(set(genelist))
genelist.remove('(Intercept)')

geneinfo = geneinfo[genelist]
geneinfo = np.log2(geneinfo)
survinfo = survinfo.join(geneinfo)
print(survinfo.columns)


def DrawCoxCruve(data, time_col, statu_col):
    plt.figure(figsize=(6, 4), dpi=300)
    cph = CoxPHFitter()
    cph.fit(survinfo, duration_col=time_col, event_col=statu_col)
    cph.print_summary()
    print(pd.concat([cph.params_, cph.confidence_intervals_], axis=1))
    print(cph.confidence_intervals_)
    cph.plot()
    return plt


plt = DrawCoxCruve(data=survinfo,
                   time_col='days_to_death',
                   statu_col='vital_status')

#plt.xlim(-1.25, 1.25)
#plt.tight_layout()
#plt.savefig('/Users/freud/Documents/MANU/BI/PANC_ECM/RESULT/ha_cox.jpg')
