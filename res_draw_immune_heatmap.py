#!/usr/bin/env python
# -*- coding: utf-8 -*-
# @Author : Ji Hongchen
# @Email : jhca.xyt@163.com
# @Last version : 2022-01-12 16:10
# @Filename : res_draw_meta_heatmap.py

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import cfg
import seaborn as sns


CFG = cfg.cfg()
db = 'TCGA'

if db == 'TCGA':
    EXP_FILE_ALL = pd.read_csv(
        CFG.datapath + 'SELECTED_MRNA.csv', index_col=[0])
    class_file = pd.read_csv('~/Documents/MANU/BI/PANC_ECM/Manuscript/Data/TYPE2_post_lasso_coef_2.csv',
                             index_col=[0])
elif db == 'CPT':
    EXP_FILE_ALL = pd.read_csv(
        CFG.datapath + 'CPTAC_MRNA_SYMBOL.csv', index_col=[0])
    EXP_FILE_ALL = pd.DataFrame(np.array(EXP_FILE_ALL).T,
                                index=EXP_FILE_ALL.columns, columns=EXP_FILE_ALL.index)
    class_file = pd.read_csv(
        '~/Documents/MANU/BI/PANC_ECM/Manuscript/Data/Clinica_CPT_TYPE2_COEF_2.csv',
        index_col=[0])[['lasso_result']]
    class_file.columns = ['lasso_class']

print(EXP_FILE_ALL)

#for i in EXP_FILE_ALL.columns:
#    exp_i_m = EXP_FILE_ALL[i].mean()
#    exp_i_s = EXP_FILE_ALL[i].std()
#    EXP_FILE_ALL[i] = (EXP_FILE_ALL[i] - exp_i_m) / exp_i_s

genefile = '~/Documents/MANU/BI/PANC_ECM/DATA/Immune_type_genes.csv'
# genefile = '~/Documents/MANU/BI/PANC_ECM/Manuscript/Data/ROS_genes.csv'
geneinfo = pd.read_csv(genefile)

pathway_columns = list(set(geneinfo['Cluster name']))
for pathway in pathway_columns:
    genelist = list(geneinfo['entrez'][geneinfo['Cluster name'] == pathway])
    if 0 in genelist:
        genelist.remove(0)
    EXP_FILE = EXP_FILE_ALL.loc[genelist]
    EXP_FILE = np.log2(EXP_FILE)
    EXP_FILE.replace(-np.inf, np.nan, inplace=True)
    EXP_FILE.fillna(EXP_FILE.mean(), inplace=True)
    high_express = list(
        class_file[class_file['lasso_class'] == 'STROMA_HIGH'].index)
    low_express = list(
        class_file[class_file['lasso_class'] == 'STROMA_LOW'].index)
    high_mat = np.array(EXP_FILE[high_express])
    low_mat = np.array(EXP_FILE[low_express])
    insert = np.ones((len(high_mat), 3)) * 0
    drawdata = np.hstack((high_mat, insert, low_mat))
    drawdata = np.hstack((high_mat, low_mat))
    for i in range(len(drawdata)):
        drawdata[i] = (drawdata[i] - drawdata[i].mean()) / drawdata[i].std()
    drawdata[:, high_mat.shape[1]: high_mat.shape[1] + 3] = 0

    figure = plt.figure(figsize=(5, 5), dpi=300)
    plt.style.use('ggplot')
    sns.set_style('whitegrid')
    sns.heatmap(drawdata,
                cbar=True,
                vmin=-1.5,
                vmax=1.5,
                cmap='bwr',
                yticklabels=genelist)
    plt.yticks(fontsize=20)
    plt.xticks([])
    plt.xlabel('')
    plt.tight_layout()
    plt.savefig('/Users/freud/Documents/MANU/BI/PANC_ECM/Manuscript/Figure/SubFigure/' +
                '_'.join(pathway.split(' ')) + '_Immune_heatmap.tiff')
