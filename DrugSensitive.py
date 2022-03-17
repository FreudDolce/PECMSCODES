#!/usr/bin/env python
# -*- coding: utf-8 -*-
# @Author : Ji Hongchen
# @Email : jhca.xyt@163.com
# @Last version : 2021-12-01 22:55
# @Filename : drugsensitive.py

import pandas as pd
import numpy as np
import os
import argparse
import cfg
from DrawBoxPlots import MergeClassFile, Platedata, DrawMutiBoxPlot

CFG = cfg.cfg()

parser = argparse.ArgumentParser()
parser.add_argument('-n', help='Number of cluster')
args = parser.parse_args()

workspace = CFG.resultpath + 'SIG_RESULT/'

mrnadata = CFG.EXPRESSION_FILE
drugs = ['Gemcitabine', 'Erlotinib', 'Docetaxel', 'Cisplatin',
         'Paclitaxel', '5-Fluorouracil', 'Thapsigargin']
if __name__ == '__main__':
    folderlist = os.listdir(workspace)
    for drug in drugs:
        print('Drug: ', drug)
        savename = drug + '-'
        items = [drug]
        for folder in folderlist:
            print('Folder: ', folder)
            nf = MergeClassFile(
                CFG.datapath + 'DrugSensitive_log2.csv',
                workspace + folder + '/ClusterResult/cluster_' +
                str(args.n) + '.csv',
                workspace + folder + '/classdict_' + str(args.n) + '.npy',
                has_inplace=True,
                has_changeid=False
            )
            pf = Platedata(nf)
            plt = DrawMutiBoxPlot(
                dataframe=pf,
                statcols=['hue', 'classname'],
                valuecol='value',
                statitem=items
            )
            plt.xticks(rotation=90)
            plt.tight_layout()
            plt.savefig(workspace + folder + '/Drug-sensitive-' +
                        savename + str(args.n) + '_box.jpg')
            plt.close()
