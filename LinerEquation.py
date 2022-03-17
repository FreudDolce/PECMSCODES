#!/usr/bin/env python
# -*- coding: utf-8 -*-
# @Author : Ji Hongchen
# @Email : jhca.xyt@163.com
# @Last version : 2021-12-20 19:37
# @Filename : LinerEquation.py

import pandas as pd
import numpy as np
import os
import cfg
from sklearn.metrics import roc_curve, auc
import matplotlib.pyplot as plt

CFG = cfg.cfg()

workspace = CFG.resultpath + 'SIG_RESULT/0.0021-str-comb_20211125/SCORE/'
FUNC = 'COL'
USE_COEF = '2'

exp_data = workspace + FUNC + '_LASSO.csv'
coef_data = workspace + FUNC + '_COEF_' + USE_COEF + '.csv'


def GetPredResult(data, coef):
    exp_info = pd.read_csv(data, index_col=[0])
    coef_info = pd.read_csv(coef_data, index_col=[0])
    coef_info.drop(coef_info[coef_info['s1'] == '.'].index, inplace=True)
    incept = float(coef_info['s1']['(Intercept)'])
    coef_info.drop(index='(Intercept)', axis=1, inplace=True)
    genelist = list(coef_info.index)
    mut_l = np.array(exp_info[genelist]).astype('float')
    genelist.append('cluster')
    exp_info = exp_info[genelist]
    mut_r = np.array(coef_info).astype('float')
    pred = (np.dot(mut_l, mut_r) + incept).reshape(-1, 1)
    exp_info['pred'] = pred
    return exp_info


def DrawRocCurve(pred, label):
    fpr, tpr, thersholds = roc_curve(pred, label)
    thershold = thersholds[np.argmax(tpr - fpr)]
    print('Thershold: ', thershold)
    roc_auc = auc(fpr, tpr)
    fig = plt.figure(figsize=(6, 6), dpi=300)
    plt.plot(fpr, tpr, 'k--',
             label='ROC (area = {0:.2f})'.format(roc_auc), lw=2)
    plt.xlim([-0.05, 1.05])  # 设置x、y轴的上下限，以免和边缘重合，更好的观察图像的整体
    plt.ylim([-0.05, 1.05])
    plt.xlabel('False Positive Rate', size=20)
    plt.ylabel('True Positive Rate', size=20)  # 可以使用中文，但需要导入一些库即字体
    plt.tick_params(labelsize=15)
    plt.legend(loc="lower right")
    plt.savefig(workspace + FUNC + '_ROC.jpg')
    return thershold


def AddLassoResult(df, pred_col, thres):
    df['lasso_class'] = 1
    df['lasso_class'][df[pred_col] >= threshold] = 2
    return df


if __name__ == '__main__':
    df = GetPredResult(exp_data, coef_data)
    threshold = DrawRocCurve(list(df['cluster']), list(df['pred']))
    ldf = AddLassoResult(df, 'pred', threshold)
    wrong_num = len(ldf[ldf['lasso_class'] != ldf['cluster']])
    print('Accuracy rate: ', (len(ldf) - wrong_num) / len(ldf))
    ldf = ldf[['cluster', 'pred', 'lasso_class']]
    ldf.to_csv(workspace + FUNC + '_POST_LASSO.csv')
