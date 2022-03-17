#!/usr/bin/env python
# -*- coding: utf-8 -*-
# @Author : Ji Hongchen
# @Email : jhca.xyt@163.com
# @Last version : 2021-10-05 10:11
# @Filename : SurvivalAnalysis.py


import pandas as pd
import numpy as np
import os
import json
import cfg

CFG = cfg.cfg()

CLIC_CSV = CFG.datapath + 'clinical_TCGA_for_R.csv'
CLAS_FOLD = CFG.CLUSTER_RESULT_PATH


def MergeClinicalAndClusterResult(clinicaldata, classificationdata):
    clinicalinfo = pd.read_csv(clinicaldata)
    clinicalinfo = MergeSurvivalTime(clinicalinfo)
    classficationinfo = pd.read_csv(CLAS_FOLD + classificationdata)
    classficationinfo['case_id'] = 'a'
    for i in classficationinfo.index:
        if 'X' in classficationinfo['Unnamed: 0'][i]:
            classficationinfo['case_id'][i] = \
                '-'.join(classficationinfo['Unnamed: 0']
                         [i].split('X')[1].split('.'))
        else:
            classficationinfo['case_id'][i] = '-'.join(
                classficationinfo['Unnamed: 0'][i].split('.'))
    classficationinfo = classficationinfo.set_index('case_id')

    clinicalinfo['class'] = 11
    for j in clinicalinfo.index:
        try:
            clinicalinfo['class'][j] = classficationinfo['x'].loc[clinicalinfo['case_id'][j]]
        except KeyError:
            pass
    clinicalinfo = clinicalinfo.drop(
        clinicalinfo[clinicalinfo['class'] == 11].index)
    return clinicalinfo


def MergeSurvivalTime(clinicalinfo):
    for i in clinicalinfo.index:
        if clinicalinfo['vital_status'][i] == 'Alive':
            clinicalinfo['days_to_death'][i] = clinicalinfo['days_to_last_follow_up'][i]
    clinicalinfo['vital_status'][clinicalinfo['vital_status'] == 'Alive'] = 1
    clinicalinfo['vital_status'][clinicalinfo['vital_status'] == 'Dead'] = 2
    return clinicalinfo


if __name__ == '__main__':
    for c in range(CFG.MIN_CLASS, CFG.MAX_CLASS + 1):
        df = MergeClinicalAndClusterResult(CLIC_CSV, 'cluster_' + str(c) + '.csv')
        df.to_csv(CLAS_FOLD + 'With_class_cluster_' + str(c) + '.csv')
