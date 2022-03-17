#!/usr/bin/env python
# -*- coding: utf-8 -*-
# @Author : Ji Hongchen
# @Email : jhca.xyt@163.com
# @Last version : 2021-10-08 20:35
# @Filename : CIBERSORT.py

import pandas as pd
import numpy as np
import os
import cfg
from scipy import stats

CFG = cfg.cfg()


def CiberPvalueTest(class_n):
    ciber_result = pd.read_csv(CFG.CIBERSORT_RESULT)
    immun_col = ciber_result.columns[1: -3]
    ciber_p_frame = pd.DataFrame(
        np.zeros((2, 22)), columns=immun_col, index=['f', 'p'])
    ciber_result['class'] = 0
    class_info = pd.read_csv(CFG.CLUSTER_RESULT_PATH + 'With_class_cluster_' + str(class_n) + '.csv',
                             index_col=[0])
    class_info = class_info.set_index('case_id')
    ciber_result['class'] = 20
    for i in ciber_result.index:
        ciber_result['class'][i] = \
            class_info['class'][ciber_result['Unnamed: 0'][i]]
    ciber_result = ciber_result.drop(
        ciber_result[ciber_result['class'] == 20].index)
    ciber_result.to_csv(CFG.CLUSTER_RESULT_PATH + 'ciber_result_class_' + str(class_n) + '.csv')
    kwargs = []
    for im in immun_col:
        for i in range(class_n):
            kwargs.append(
                list(ciber_result[ciber_result['class'] == i + 1][im]))
        f, p = stats.f_oneway(*kwargs)
        ciber_p_frame[im]['f'] = f
        ciber_p_frame[im]['p'] = p
    for c in range(1, class_n + 1):
        ciber_p_frame.loc[c] = 0
        for im in immun_col:
            ciber_p_frame[im][c] = \
                    ciber_result[ciber_result['class'] == c][im].quantile()
    ciber_p_frame.to_csv(CFG.CLUSTER_RESULT_PATH +
                         'ciber_result_p_class_' + str(class_n) + '.csv')


if __name__ == '__main__':
    for n in range(CFG.MIN_CLASS, CFG.MAX_CLASS + 1):
        CiberPvalueTest(n)
