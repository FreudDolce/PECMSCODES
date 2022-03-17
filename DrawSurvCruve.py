#!/usr/bin/env python
# -*- coding: utf-8 -*-
# @Author : Ji Hongchen
# @Email : jhca.xyt@163.com
# @Last version : 2021-10-20 22:23
# @Filename : DrawSurvCruve.py

import pandas as pd
from lifelines import KaplanMeierFitter as KMF
from lifelines.statistics import logrank_test
import os
import matplotlib.pyplot as plt
import cfg

CFG = cfg.cfg()

DRAW_RESUTL = CFG.resultpath + 'SIG_RESULT/'

def DrawSurvivalplot(dataframe, has_plot=True, class_col = 'class'):
    kmf = KMF()
    fig = plt.figure(figsize=(4.5, 4), dpi=300)
    classes = list(set(dataframe[class_col]))
    n_classes = len(classes)
    median_surv = {}
    p_frame = pd.DataFrame(columns=classes, index=classes)
    for c in classes:
        statframe = dataframe[dataframe[class_col] == c]
        _T_1 = statframe['days_to_death'].astype('int')
        _E_1 = statframe['vital_status'].astype('int')
        kmf.fit(_T_1, event_observed=_E_1,
                label='Class-' + str(c) + ' (n = ' + str(len(statframe)) + ').')
        median_surv[str(c)] = kmf.median_survival_time_
        median_surv_list = sorted(
                median_surv.items(), key=lambda d: d[1], reverse=True)
        for sub_c in classes:
            if sub_c != c:
                sub_statframe = dataframe[dataframe[class_col] == sub_c]
                _T_2= sub_statframe['days_to_death'].astype('int')
                _E_2 = sub_statframe['vital_status'].astype('int')
                lr = logrank_test(_T_1, _T_2, _E_1, _E_2, alpha=0.95)
                p_frame[c].loc[sub_c] = lr.p_value
    if has_plot == True:
        classes = list(set(dataframe[class_col]))
        for c in classes:
            statframe = dataframe[dataframe[class_col] == c]
            _T_1 = statframe['days_to_death']
            _E_1 = statframe['vital_status']
            kmf.fit(_T_1, event_observed=_E_1,
                    label=str(c) + ' (n = ' + str(len(statframe)) + ').')
            kmf.plot(ci_alpha=0)
        return (p_frame, plt)

if __name__ == '__main__':
    result_dirs = os.listdir(DRAW_RESUTL)
    for result_dir in result_dirs:
        result_path = DRAW_RESUTL + result_dir + '/ClusterResult/'
        for result in os.listdir(result_path):
            if (('.csv' in result) and ('With_class_cluster' in result)):
                df = pd.read_csv(result_path + result)
                survframe, plt = DrawSurvivalplot(df)
                #survframe.to_csv(result_path + result, index=False)
                print(survframe)
                plt.savefig(result_path + result.split('.')[0] + '.jpg')
        print(result_path + ' ' + result_path + ' saved.')
        plt.close()
