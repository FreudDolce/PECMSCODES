#!/usr/bin/env python
# -*- coding: utf-8 -*-
# @Author : Ji Hongchen
# @Email : jhca.xyt@163.com
# @Last version : 2021-12-22 18:01
# @Filename : pl_surv_analysis.py

import numpy as np
import pandas as pd
import cfg
from DrawSurvCruve import DrawSurvivalplot
from SurvivalAnalysis import MergeSurvivalTime
import os

CFG = cfg.cfg()

clinicalinfo = pd.read_csv(CFG.datapath + 'CPT_clinical.csv', index_col=[0])
classinfo = pd.read_csv(CFG.lassoresultpath + 'LASSO_CPT_COMB.csv', index_col=[0])


def GetClinicalClass(clinicaldata, classdata, class_col):
    #clinicalinfo = MergeSurvivalTime(clinicaldata)
    print(clinicalinfo)
    clinicalinfo['pl_class'] = 'a'
    for i in clinicalinfo.index:
        try:
            clinicalinfo['pl_class'].loc[i] = classdata[class_col][i]
        except KeyError:
            pass
    clinicalinfo.drop(clinicalinfo[clinicalinfo['pl_class'] == 'a'].index,
                      inplace=True)
    return clinicalinfo


if __name__ == '__main__':
    class_by = 'COMB'
    cf = GetClinicalClass(clinicalinfo, classinfo, class_by)
    sf, plt = DrawSurvivalplot(cf, has_plot=True, class_col='pl_class')
    plt.savefig(CFG.lassoresultpath + 'surv_cpt' + class_by + '.jpg')
    sf.to_csv(CFG.lassoresultpath + 'surv_cpt' + class_by + '.csv')
