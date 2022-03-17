#!/usr/bin/env python
# -*- coding: utf-8 -*-
# @Author : Ji Hongchen
# @Email : jhca.xyt@163.com
# @Last version : 2022-01-20 20:45
# @Filename : res_draw_corr_heat_map.py

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import os

path = '/Users/freud/Documents/MANU/BI/PANC_ECM/Manuscript/Data/'
folders = os.listdir(path)

for f in folders:
    if 'Draw_corr_' in f:
        print(f)
        drawframe = pd.read_csv(path + f, index_col=[0])
        X = list(drawframe.columns)
        Y = list(drawframe.index)

        figure = plt.figure(figsize=(10, 5), dpi=300)
        plt.style.use('ggplot')
        sns.set_style('whitegrid')

        sns.heatmap(drawframe,
                    cbar=True,
                    vmin=-0.6,
                    vmax=0.6,
                    cmap='bwr',
                    xticklabels=X,
                    yticklabels=Y)
        ax = plt.gca()
        bwith = 3  # 边框宽度设置为2
        ax = plt.gca()  # 获取边框
        ax.spines['bottom'].set_linewidth(bwith)
        ax.spines['left'].set_linewidth(bwith)
        ax.spines['top'].set_linewidth(bwith)
        ax.spines['right'].set_linewidth(bwith)
        plt.tick_params(width=2.5, length=5.5)
        plt.yticks(fontproperties='Arial', size=25)
        plt.xticks(fontproperties='Arial', size=25, rotation=90)
        plt.legend(bbox_to_anchor=(1, 1.4))
        plt.gcf().subplots_adjust(left=0.3, right=0.92, bottom=0.4, top=0.92)
        plt.savefig('/Users/freud/Documents/MANU/BI/PANC_ECM/Manuscript/Figure/SubFigure/' +
                    f.split('.')[0] + '.tiff')
