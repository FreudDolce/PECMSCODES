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

SIGFILE = 'TYPE2'
CFG = cfg.cfg()

EXP_FILE = '~/Documents/MANU/BI/PANC_ECM/DATA/SELECTED_MRNA_SYMBOL.csv'
SURV_FILE = '~/Documents/MANU/BI/PANC_ECM/Manuscript/Data/TCGA_ClinicalwithCluster.csv'
expinfo = pd.read_csv(EXP_FILE, index_col=[0])
geneinfo = pd.DataFrame(np.array(expinfo).T,
                        index=expinfo.columns,
                        columns=expinfo.index)
survinfo = pd.read_csv(SURV_FILE)
survinfo = survinfo[['case_id', 'days_to_death', 'vital_status']]
survinfo.set_index('case_id', drop=True, inplace=True)

# genes = pd.read_csv('~/Documents/MANU/BI/PANC_ECM/Manuscript/Data/' + SIGFILE + '_limma_sig.csv',
#                    index_col=[0])
#genelist = list(genes.index)
genelist = list(pd.read_csv('~/Documents/MANU/BI/PANC_ECM/Manuscript/Data/TYPE2_lasso_coef2.csv',
                         index_col=[0]).index)
genelist.remove('(Intercept)')
genelist = list(set(genelist))

# genelist = list(pd.read_csv('~/Documents/MANU/BI/PANC_ECM/Manuscript/Data/TYPE2_lasso_coef1.csv',
#                            index_col=[0]).index)
# genelist.remove('(Intercept)')
geneinfo = geneinfo[genelist]
print(geneinfo)
geneinfo = np.log2(geneinfo)
survinfo = survinfo.join(geneinfo)


def SPCox(data, time_col, statu_col):
    for gene in genelist:
        cols = ['days_to_death', 'vital_status']
        cols.append(gene)
        singleinfo = data[cols]
        cph = CoxPHFitter()
        cph.fit(singleinfo, duration_col=time_col, event_col=statu_col)
        item_frame = pd.concat(
            [cph.params_, cph.confidence_intervals_], axis=1)
        try:
            totalframe = pd.concat([totalframe, item_frame], axis=0)
        except NameError:
            totalframe = item_frame
    totalframe.columns = ['coef', 'lower', 'upper']
    totalframe['sig'] = totalframe['upper'] * totalframe['lower']
    totalframe = totalframe[totalframe['sig'] >= 0]
    return totalframe


def DrawCoxCruve(data, time_col, statu_col):
    plt.figure(figsize=(5, 10), dpi=300)
    cph = CoxPHFitter()
    cph.fit(data, duration_col=time_col, event_col=statu_col)
    output_frame = pd.concat([cph.params_, cph.confidence_intervals_], axis=1)
    output_frame.columns = ['coef', 'lower', 'upper']
    output_frame['sig'] = output_frame['upper'] * output_frame['lower']
    output_frame.to_csv(
        '~/Documents/MANU/BI/PANC_ECM/Manuscript/Data/' + SIGFILE + '_Part_MutiCox.csv')
    print(output_frame)
    output_frame = output_frame[output_frame['sig'] >= 0]
    cph.plot()
    return (plt, output_frame)


if __name__ == '__main__':
    sigframe = SPCox(data=survinfo,
                     time_col='days_to_death',
                     statu_col='vital_status')
    sigframe.to_csv(
        '~/Documents/MANU/BI/PANC_ECM/Manuscript/Data/Part_SingleCoxSig.csv')
    siglist = list(sigframe.index)
    siglist.extend(['days_to_death', 'vital_status'])
    mutiinfo = survinfo[siglist]
    plt, output_frame = DrawCoxCruve(data=mutiinfo,
                                     time_col='days_to_death',
                                     statu_col='vital_status')
    print(output_frame)
    output_frame = output_frame[['coef']]
    output_frame.columns = ['s1']
    #output_frame.to_csv(
    #    '~/Documents/MANU/BI/PANC_ECM/Manuscript/Data/' + SIGFILE + '_lasso_coef3.csv')
    plt.tight_layout()
    plt.savefig(
        '/Users/freud/Documents/MANU/BI/PANC_ECM/Manuscript/Figure/SubFigure/' + SIGFILE + '_Part_Muticox.jpg')
