#!/usr/bin/env python
# -*- coding: utf-8 -*-
# @Author : Ji Hongchen
# @Email : jhca.xyt@163.com
# @Last version : 2022-03-10 10:29
# @Filename : res_draw_ihc_count.py

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from scipy import stats

Drawlist = ['CD8', 'CD31', 'PD-L1', 'GLUT-1']
EXP = pd.read_csv('~/Documents/MANU/BI/PANC_ECM/Manuscript/Data/IHC_count.csv',
                  index_col=[0])

for gene in Drawlist:
    gene_in_high = np.array(EXP[gene][EXP['PECMS_class'] == 1])
    gene_in_low = np.array(EXP[gene][EXP['PECMS_class'] == 0])
    p = stats.mannwhitneyu(gene_in_high, gene_in_low, use_continuity=True)
    print(gene, ':', p.pvalue)

    figure = plt.figure(figsize=(2.35, 4), dpi=300)
    sns.swarmplot(x='PECMS_class',
                  y=gene,
                  data=EXP,
                  palette='deep',
                  alpha=0.8,
                  size=3,
                  order=[1, 0]
                  )
    sns.boxplot(x='PECMS_class',
                y=gene,
                data=EXP,
                color='white',
                linewidth=2.5,
                width=0.6,
                saturation=1,
                order=[1, 0]
                )
    plt.title('P = ' + str(round(p.pvalue, 3)))
    bwith = 2.5  # 边框宽度设置为2
    ax = plt.gca()  # 获取边框
    ax.spines['bottom'].set_linewidth(bwith)
    ax.spines['left'].set_linewidth(bwith)
    ax.spines['top'].set_linewidth(bwith)
    ax.spines['right'].set_linewidth(bwith)
    plt.tick_params(width=2.5, length=5.5)
    plt.xlabel(' ')
    plt.ylabel(' ')
    plt.xticks([])
    plt.yticks(fontproperties='Arial', size=25)
    plt.legend(bbox_to_anchor=(1, 1.4))
    plt.gcf().subplots_adjust(left=0.4, right=0.92)
    # plt.tight_layout()
    plt.savefig(
        '/Users/freud/Documents/MANU/BI/PANC_ECM/Manuscript/Figure/SubFigure/IHC_exp_' +
        str(gene) + '_box.tiff', transparent=True)
