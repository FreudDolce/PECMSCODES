#!/usr/bin/env python
# -*- coding: utf-8 -*-
# @Author : Ji Hongchen
# @Email : jhca.xyt@163.com
# @Last version : 2021-11-20 21:35
# @Filename : DrawHeatMap.py

import pandas as pd
import numpy as np
import os
import matplotlib.pyplot as plt
import seaborn as sns
import cfg
from DrawBoxPlots import _ExchangePatientID

CFG = cfg.cfg()
GSVA_BLANK_LINE = 2
MRNA_BLANK_LINE = 10


def DrawHeatMap(dataframe, l_width, has_xlabel=True):
    plt.figure(figsize=(6, 6), dpi=300)
    sns.heatmap(dataframe,
                #vmax=show_max,
                #vmin=show_min,
                linewidths=l_width,
                cmap='RdBu')
    return (plt)


def GetClassDict(classframe):
    cdit = {}
    for c in range(len(list(set(classframe['x'])))):
        cdit[c + 1] = list(classframe[classframe['x'] == c + 1]['Unnamed: 0'])
    ncdict = {}
    for n in cdit:
        ncdict[n] = []
        for i in cdit[n]:
            ncdict[n].append(_ExchangePatientID(i))
    return ncdict


if __name__ == '__main__':
    gsvaresult = pd.read_csv(
        CFG.datapath + 'drug_sensitive_gsva_result.csv', index_col=[0])
    gsva_color = [gsvaresult.quantile(0.1).median(),
                  gsvaresult.quantile(0.3).median(),
                  gsvaresult.quantile(0.6).median(),
                  gsvaresult.quantile(0.9).median()]
    n_cluster = [3, 4]
    folderlist = os.listdir(CFG.resultpath + 'SIG_RESULT/')
    for folder in folderlist:
        print('======================================================')
        print('Folder: ', folder)
        for n in n_cluster:
            print('------------------------------------------------------')
            print('Number of classes: ', n)
            classframe = pd.read_csv(
                CFG.resultpath + 'SIG_RESULT/' + folder + '/ClusterResult/cluster_' + str(n) + '.csv')
            classdict = GetClassDict(classframe)
            difindict = np.load(
                CFG.resultpath + 'SIG_RESULT/' + folder +
                '/classdict_' + str(n) + '.npy',
                allow_pickle=True).item()
            gsva_combined_cols = []
            cols = []
            for c in classdict:
                gsva_combined_cols.append(difindict[c])
                cols.extend(classdict[c])
            gsvaresult_draw = pd.DataFrame(
                columns=gsva_combined_cols, index=gsvaresult.index)
            for c in classdict:
                for ind in gsvaresult_draw.index:
                    gsvaresult_draw[difindict[c]].loc[ind] = \
                        gsvaresult[classdict[c]].loc[ind].median()
            #for G_BLANK in range(GSVA_BLANK_LINE):
            #    gsvaresult_draw.loc[len(gsvaresult_draw)] = 0.0
            #for g_line in range(len(gsvaresult_draw) - int(GSVA_BLANK_LINE / 2), len(gsvaresult_draw)):
            #    for c in classdict:
            #        gsvaresult_draw.loc[g_line][difindict[c]
            #                                    ] = gsva_color[c - 1]
            gsvaresult_draw = gsvaresult_draw.astype('float64')
            print(gsvaresult_draw)
            plt = DrawHeatMap(gsvaresult_draw,
                              #show_max=0.4,
                              #show_min=-0.4,
                              l_width=0.1)
            plt.savefig(CFG.resultpath + 'SIG_RESULT/' + folder + '/ClusterResult/drug_sens_gsva_m_' + str(n) + '.jpg',
                        bbox_inches='tight')
            plt.close()

    gsvaresult = pd.read_csv(
        CFG.datapath + 'func_gsva_result.csv', index_col=[0])
    gsva_color = [gsvaresult.quantile(0.1).median(),
                  gsvaresult.quantile(0.3).median(),
                  gsvaresult.quantile(0.6).median(),
                  gsvaresult.quantile(0.9).median()]
    n_cluster = [3, 4]
    folderlist = os.listdir(CFG.resultpath + 'SIG_RESULT/')
    for folder in folderlist:
        print('======================================================')
        print('Folder: ', folder)
        for n in n_cluster:
            print('------------------------------------------------------')
            print('Number of classes: ', n)
            classframe = pd.read_csv(
                CFG.resultpath + 'SIG_RESULT/' + folder + '/ClusterResult/cluster_' + str(n) + '.csv')
            classdict = GetClassDict(classframe)
            difindict = np.load(
                CFG.resultpath + 'SIG_RESULT/' + folder +
                '/classdict_' + str(n) + '.npy',
                allow_pickle=True).item()
            gsva_combined_cols = []
            cols = []
            for c in classdict:
                gsva_combined_cols.append(difindict[c])
                cols.extend(classdict[c])
            gsvaresult_draw = pd.DataFrame(
                columns=gsva_combined_cols, index=gsvaresult.index)
            for c in classdict:
                for ind in gsvaresult_draw.index:
                    gsvaresult_draw[difindict[c]].loc[ind] = \
                        gsvaresult[classdict[c]].loc[ind].median()
            #for G_BLANK in range(GSVA_BLANK_LINE):
            #    gsvaresult_draw.loc[len(gsvaresult_draw)] = 0.0
            #for g_line in range(len(gsvaresult_draw) - int(GSVA_BLANK_LINE / 2), len(gsvaresult_draw)):
            #    for c in classdict:
            #        gsvaresult_draw.loc[g_line][difindict[c]
            #                                    ] = gsva_color[c - 1]
            gsvaresult_draw = gsvaresult_draw.astype('float64')
            print(gsvaresult_draw)
            plt = DrawHeatMap(gsvaresult_draw,
                              #show_max=0.4,
                              #show_min=-0.4,
                              l_width=0.1)
            plt.savefig(CFG.resultpath + 'SIG_RESULT/' + folder + '/ClusterResult/test_func_sens_gsva_m_' + str(n) + '.jpg',
                        bbox_inches='tight')
            plt.close()

            """
            rnaexp_draw = rnaexp[cols]
            for R_BLANK in range(MRNA_BLANK_LINE):
                rnaexp_draw.loc[len(rnaexp_draw)] = 1
            for r_line in range(len(rnaexp_draw) - int(MRNA_BLANK_LINE / 2), len(rnaexp_draw)):
                for c in classdict:
                    rnaexp_draw.loc[r_line][classdict[c]] = rnaexp_color[c - 1]
            rnaexp_draw = np.log2(rnaexp_draw)
            plt = DrawHeatMap(rnaexp_draw,
                              show_max=0,
                              show_min=-12,
                              l_width=0)
            plt.savefig(CFG.resultpath + 'SIG_RESULT/' + folder + '/ClusterResult/drug_sens_mrna_' + str(n) + '.jpg',
                        bbox_inches='tight')
            plt.close()

            # for blank in range(BLANK_LINE):
            #    rnaexp.loc[len(rnaexp) + 1] = 0
            # for clasline in range(len())
            """
