#!/usr/bin/env python
# -*- coding: utf-8 -*-
# @Author : Ji Hongchen
# @Email : jhca.xyt@163.com
# @Last version : 2022-01-06 00:07
# @Filename : res_draw_tmb.py
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from scipy import stats

COEF = '2'
item = 'TYPE2'
db = 'CPT'

tmb_data = pd.read_csv('~/Documents/MANU/BI/PANC_ECM/Manuscript/Data/TMB_' + db + '.csv',
                       index_col=[0])
tmb_data = tmb_data[['Tumor_Sample_Barcode', 'total_perMB']]
tmb_data['total_perMB'] = np.log10(tmb_data['total_perMB'])
tmb_data.replace(-np.inf, np.nan, inplace=True)
tmb_data.dropna(how='any', inplace=True)
clincial_data = pd.read_csv('~/Documents/MANU/BI/PANC_ECM/Manuscript/Data/Clinica_' +
                            db + '_' + item + '_COEF_' + COEF + '.csv',
                            index_col=[0])
clincial_data = pd.merge(clincial_data, tmb_data, on='Tumor_Sample_Barcode')
print(clincial_data)


def DrawBoxMap(class_by, value, dataframe, lw, draw_order=[]):
    high_value = dataframe[value][dataframe[class_by] == 'STROMA_HIGH']
    low_value = dataframe[value][dataframe[class_by] == 'STROMA_LOW']
    p = stats.mannwhitneyu(high_value, low_value)
    figure = plt.figure(figsize=(4.3, 8), dpi=300)
    sns.boxplot(x=class_by, y=value,
                data=dataframe, color='white',
                linewidth=lw,
                width=0.7, palette='deep',
                order=draw_order)
    plt.xticks(rotation='vertical', size=20)
    plt.title('P: ' + str(round(p.pvalue, 5)))
    plt.yticks(size=20)
    return plt


def DrawVioMap(class_by, value, dataframe, lw, draw_order=[]):
    high_value = dataframe[value][dataframe[class_by] == 'STROMA_HIGH']
    print(high_value)
    low_value = dataframe[value][dataframe[class_by] == 'STROMA_LOW']
    p = stats.mannwhitneyu(high_value, low_value)
    print(p.pvalue)
    figure = plt.figure(figsize=(4.3, 8), dpi=300)
    sns.violinplot(x=class_by, y=value,
                   data=dataframe, color='white',
                   linewidth=lw,
                   width=0.7, palette='deep',
                   order=draw_order)
    plt.xticks(rotation='vertical', size=20)
    plt.yticks(size=20)
    plt.title('P: ' + str(round(p.pvalue, 5)))
    return plt


plt = DrawVioMap(class_by='lasso_result',
                 value='total_perMB',
                 dataframe=clincial_data,
                 lw=2.5,
                 draw_order=['STROMA_HIGH', 'STROMA_LOW'])
# sns.swarmplot(x='lasso_result', y='total_perMB',
#              data=clincial_data, color='black',
#              order=['STROMA-HIGH', 'STROMA-LOW'],
#              size=7)
bwith = 4 #边框宽度设置为2
ax = plt.gca()#获取边框
ax.spines['bottom'].set_linewidth(bwith)
ax.spines['left'].set_linewidth(bwith)
ax.spines['top'].set_linewidth(bwith)
ax.spines['right'].set_linewidth(bwith)
plt.tick_params(width=4, length=5.5)
plt.ylim(-3, 3)
plt.gcf().subplots_adjust(left=0.22, right=0.92)
plt.savefig('/Users/freud/Documents/MANU/BI/PANC_ECM/Manuscript/Figure/SubFigure/TMB_cluster_box_' +
            db + '_' + item + '_COEF_' + COEF + '.tiff', transparent=True)
plt.close()
