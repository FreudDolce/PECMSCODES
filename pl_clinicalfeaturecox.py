#!/usr/bin/env python
# -*- coding: utf-8 -*-
# @Author : Ji Hongchen
# @Email : jhca.xyt@163.com
# @Last version : 2021-12-28 12:14
# @Filename : pl_coxanalysis.py

import pandas as pd
import numpy as np
from lifelines import CoxPHFitter
import cfg
import matplotlib.pyplot as plt
CFG = cfg.cfg()

SURV_FILE = CFG.lassoresultpath + 'CLINICAL_WITH_CLASS.csv'
survinfo = pd.read_csv(SURV_FILE)
survinfo = survinfo[survinfo.columns[2:]]
survinfo['comb_pred'] = survinfo['col_pred'] + survinfo['ha_pred']
del survinfo['COMB']
del survinfo['COL']
del survinfo['HA']
del survinfo['race']
del survinfo['col_pred']
del survinfo['ha_pred']
#del survinfo['gender']
print(survinfo.columns)


def DrawCoxCruve(data, time_col, statu_col):
    plt.figure(figsize=(6, 4), dpi=300)
    cph = CoxPHFitter()
    cph.fit(survinfo, duration_col=time_col, event_col=statu_col)
    cph.print_summary()
    cph.data
    cph.plot()
    return plt


if __name__ == '__main__':
    DrawCoxCruve(data=survinfo,
                 time_col='days_to_death',
                 statu_col='vital_status')
    #plt.xlim(-1.25, 1.25)
    plt.tight_layout()
    # plt.show()
    plt.savefig(CFG.lassoresultpath + 'clinic_cox.jpg')
