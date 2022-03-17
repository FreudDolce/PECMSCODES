#!/usr/bin/env python
# -*- coding: utf-8 -*-
# @Author : Ji Hongchen
# @Email : jhca.xyt@163.com
# @Last version : 2022-01-08 11:14
# @Filename : res_getdataframeforcox.py

import pandas as pd
import numpy as np
import os
import cfg
import matplotlib.pyplot as plt
from lifelines import CoxPHFitter

filelist = os.listdir(
    '/Users/freud/Documents/MANU/BI/PANC_ECM/Manuscript/Data/')
REP = {'ajcc_pathologic_stage': {'I': 0, 'IA': 0, 'IB': 0, 'IIA': 0, 'IIB': 0,
                                 'III': 1, 'IIIB': 1, 'IV': 1}}
TYPE = 'TYPE2'
coef = '2'
db = 'COMB'

surv_file = pd.read_csv('~/Documents/MANU/BI/PANC_ECM/Manuscript/Data/Clinica_' + db +
                        '_' + TYPE + '_COEF_' + coef + '.csv',
                        index_col=[0])
surv_file = surv_file[surv_file['ajcc_pathologic_stage'].isin(
    REP['ajcc_pathologic_stage'])]
for item in REP['ajcc_pathologic_stage']:
    surv_file['ajcc_pathologic_stage'][surv_file['ajcc_pathologic_stage'] == item] = \
        REP['ajcc_pathologic_stage'][item]
print('Len: ', len(surv_file[surv_file['tissue_or_organ_of_origin'] == 1]))

surv_file['lasso_result'][surv_file['lasso_result'] == 'STROMA-HIGH'] = 1
surv_file['lasso_result'][surv_file['lasso_result'] == 'STROMA_HIGH'] = 1
surv_file['lasso_result'][surv_file['lasso_result'] == 'STROMA-LOW'] = 0
surv_file['lasso_result'][surv_file['lasso_result'] == 'STROMA_LOW'] = 0
print(surv_file)

surv_item = ['days_to_birth', 'ajcc_pathologic_stage', 'gender',
             'tissue_or_organ_of_origin', 'alcohol_history',
             'cigarettes_per_day', 'lasso_result']


def SPCox(data, time_col, statu_col):
    for item in surv_item:
        cols = ['days_to_death', 'vital_status']
        cols.append(item)
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
    totalframe.to_csv('~/Documents/MANU/BI/PANC_ECM/Manuscript/Data/' +
                      TYPE + '_' + db + '_clinicla_SingleCox.csv')
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
        '~/Documents/MANU/BI/PANC_ECM/Manuscript/Data/' + TYPE + '_' + db + '_clinical_MutiCox.csv')
    output_frame = output_frame[output_frame['sig'] >= 0]
    cph.plot()
    return (plt, output_frame)


if __name__ == '__main__':
    sigframe = SPCox(data=surv_file,
                     time_col='days_to_death',
                     statu_col='vital_status')
    sigframe.to_csv(
        '~/Documents/MANU/BI/PANC_ECM/Manuscript/Data/' + TYPE + '_' + db + '_clinicla_SingleCoxSig.csv')
    print(sigframe)
    siglist = list(sigframe.index)
    siglist.extend(['days_to_death', 'vital_status'])
    mutiinfo = surv_file[siglist]
    print(mutiinfo)
    plt, output_frame = DrawCoxCruve(data=mutiinfo,
                                     time_col='days_to_death',
                                     statu_col='vital_status')
    #output_frame = output_frame[['coef']]
    #output_frame.columns = ['s1']
    #output_frame.loc['(Intercept)'] = 0
    # output_frame.to_csv(
    #    '~/Documents/MANU/BI/PANC_ECM/Manuscript/Data/' + SIGFILE + '_lasso_coef3.csv')
    print(output_frame)
    plt.close()
