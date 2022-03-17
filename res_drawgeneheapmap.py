#!/usr/bin/env python
# -*- coding: utf-8 -*-
# @Author : Ji Hongchen
# @Email : jhca.xyt@163.com
# @Last version : 2022-01-10 17:57
# @Filename : res_drawgeneheapmap.py

import pandas as pd
import numpy as np
import cfg
import seaborn as sns
import matplotlib.pyplot as plt

CFG = cfg.cfg()
item = 'TYPE2'
db = 'TCGA'

if db == 'TCGA':
    EXP_FILE = pd.read_csv(
        CFG.datapath + 'SELECTED_MRNA_SYMBOL.csv', index_col=[0])
    class_file = pd.read_csv('~/Documents/MANU/BI/PANC_ECM/Manuscript/Data/TYPE2_post_lasso_coef_2.csv',
                             index_col=[0])
elif db == 'CPT':
    EXP_FILE = pd.read_csv(
        CFG.datapath + 'CPTAC_MRNA_SYMBOL.csv', index_col=[0])
    EXP_FILE = pd.DataFrame(np.array(EXP_FILE).T,
                            index=EXP_FILE.columns, columns=EXP_FILE.index)
    class_file = pd.read_csv(
        '~/Documents/MANU/BI/PANC_ECM/Manuscript/Data/Clinica_CPT_TYPE2_COEF_2.csv', 
        index_col=[0])[['lasso_result']]
    class_file.columns = ['lasso_class']
genelist = list(pd.read_csv('~/Documents/MANU/BI/PANC_ECM/Manuscript/Data/TYPE2_lasso_coef2.csv',
                            index_col=[0]).index)
genelist.remove('(Intercept)')
EXP_FILE = EXP_FILE.loc[genelist]
EXP_FILE = np.log2(EXP_FILE)
high_express = list(
    class_file[class_file['lasso_class'] == 'STROMA_HIGH'].index)
low_express = list(class_file[class_file['lasso_class'] == 'STROMA_LOW'].index)

high_mat = np.array(EXP_FILE[high_express])
low_mat = np.array(EXP_FILE[low_express])

insert = np.ones((len(high_mat), 5)) * EXP_FILE.max().max()

drawdata = np.hstack((high_mat, insert, low_mat))
for i in range(len(drawdata)):
    drawdata[i] = (drawdata[i] - drawdata[i].mean()) / drawdata[i].std()

drawdata[:, high_mat.shape[1]: high_mat.shape[1] + 5] = 0


figure = plt.figure(figsize=(9, 2.7), dpi=300)
plt.style.use('ggplot')
sns.set_style('whitegrid')
sns.heatmap(drawdata,
            cbar=True,
            vmin=-2.01,
            vmax=2.01,
            cmap='bwr',
            yticklabels=genelist)
plt.yticks(fontsize=20)
plt.xticks([])
plt.xlabel('')
plt.tight_layout()
plt.savefig(
    '/Users/freud/Documents/MANU/BI/PANC_ECM/Manuscript/Figure/SubFigure/Geneexpress_' + db + '_' + item + '.tiff')
