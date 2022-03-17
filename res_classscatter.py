#!/usr/bin/env python
# -*- coding: utf-8 -*-
# @Author : Ji Hongchen
# @Email : jhca.xyt@163.com
# @Last version : 2021-12-21 22:23
# @Filename : pl_classscater.py

import pandas as pd
import cfg
import matplotlib.pyplot as plt
import numpy as np

CFG = cfg.cfg()


def sigmoid(x):
    sig = 1 / (1 + np.exp(-x))
    return sig


SIG = False
COLOR_POOL = ['crimson', 'darkorange', 'darkgreen', 'royalblue', 'teal']
COL_THERS = 0.5414213634821112
HA_THERS = 0.4661228855009054

if SIG == True:
    COL_THERS = sigmoid(COL_THERS)
    HA_THERS = sigmoid(HA_THERS)


def DrawDifScatter(df, axis_x, axis_y, hue, color, size):
    fig = plt.figure(figsize=(6, 6), dpi=300)
    classlist = list(set(df[hue]))
    if SIG == True:
        df[axis_x] = sigmoid(df[axis_x])
        df[axis_y] = sigmoid(df[axis_y])
    plt.axhline(HA_THERS, ls=':', c='gray', lw=1)
    plt.axvline(COL_THERS, ls=':', c='gray', lw=1)
    for i in range(len(classlist)):
        drawdf = df[df[hue] == classlist[i]]
        plt.scatter(drawdf[axis_x], drawdf[axis_y],
                    c=color[i], s=size, alpha=0.5)
    plt.xlim([-1.2 + COL_THERS, 1.2 + COL_THERS])
    plt.ylim([-1.2 + HA_THERS, 1.2+ HA_THERS])
    plt.xlabel(axis_x)
    plt.ylabel(axis_y)
    return plt


if __name__ == '__main__':
    tcga_data = pd.read_csv('~/Documents/MANU/BI/PANC_ECM/Manuscript/Data/LASSO_CLASS_COMB_COEF_1.csv')
    plt = DrawDifScatter(df=tcga_data,
                         axis_x='col_pred',
                         axis_y='ha_pred',
                         hue='ori_class',
                         color=COLOR_POOL,
                         size=20)
    plt.savefig('/Users/freud/Documents/MANU/BI/PANC_ECM/Manuscript/Figures/tcga_class_scatter_coef_1.jpg')
