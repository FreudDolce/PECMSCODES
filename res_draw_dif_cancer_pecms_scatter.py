#!/usr/bin/env python
# -*- coding: utf-8 -*-
# @Author : Ji Hongchen
# @Email : jhca.xyt@163.com
# @Last version : 2022-02-22 23:38
# @Filename : res_draw_dif_cancer_pecms_scatter.py

import pandas as pd
import numpy as np
import os
import matplotlib.pyplot as plt

path = '/Users/freud/Documents/MANU/BI/PANC_ECM/DATA/TCGA-mRNA/'

gene = ['COL17A1',
        'AREG',
        'KLHL32',
        'CDA',
        'POSTN',
        'SLC2A1',
        'FN1',
        'INHBA']
coef = [0.0100817579921901,
        0.0220009100877722,
        -0.00790381705330212,
        0.0091173826150082,
        0.0167918929935681,
        0.0408528252763084,
        0.0053397646458669700,
        0.0131436954681402, ]
intercept = -1.75536346491703

ori_mrna = pd.read_csv('~/Documents/MANU/BI/PANC_ECM/DATA/SELECTED_MRNA_SYMBOL.csv',
                       index_col=[0])
ori_mrna = pd.DataFrame(np.array(ori_mrna).T,
                        index=ori_mrna.columns,
                        columns=ori_mrna.index)
ori_mrna = np.log2(ori_mrna)
ori_mrna_m = np.array(ori_mrna).mean()
ori_mrna_s = np.array(ori_mrna).std()
ori_mrna = ori_mrna[gene]


filelist = os.listdir(path)
cancers = []
boxes = []
boxes_m = []

for f in range(len(filelist)):
    filename = filelist[f].split('.')[1]
    exp_info = pd.read_csv(path + filelist[f], sep='\t', index_col=[0])
    exp_info = pd.DataFrame(np.array(exp_info).T,
                            index=exp_info.columns,
                            columns=exp_info.index)
    exp_info_m = np.array(exp_info).mean()
    exp_info_s = np.array(exp_info).std()
    exp_info = (exp_info - exp_info_m) / exp_info_s
    exp_info = exp_info * ori_mrna_s + ori_mrna_m
    exp_info = exp_info[gene]

    pecms = np.array(exp_info) @ coef + intercept
    cancers.append(filename)
    boxes.append(pecms)
    boxes_m.append(pecms.mean())
boxes_m, boxes, cancers = zip(*sorted(zip(boxes_m, boxes, cancers)))
figure = plt.figure(figsize=(3, 6), dpi=250)
plt.boxplot(boxes, labels=cancers, vert=False)
plt.xticks(rotation=90)
plt.gcf().subplots_adjust(left=0.3, right=0.95)
plt.savefig(
    '/Users/freud/Documents/MANU/BI/PANC_ECM/Manuscript/Figure/SubFigure/pecms_scatter_n.tiff',
    transparent=True)
