#!/usr/bin/env python
# -*- coding: utf-8 -*-
# @Author : Ji Hongchen
# @Email : jhca.xyt@163.com
# @Last version : 2021-12-22 21:59
# @Filename : pl_boxplot.py

import pandas as pd
import numpy as np
import cfg
from DrawBoxPlots import DrawMutiBoxPlot, Platedata, DrawMutiVioPlot, MergeClassFile, _ExchangePatientID


CFG = cfg.cfg()

#DATA_FILE = '/Users/freud/Documents/data/PANC_TCGA/DATA/estimate_score_for_plot.csv'
DATA_FILE = '/Users/freud/Documents/data/PANC_TCGA/DATA/CIBERSORT_RESULT_for_plot.csv'
#DATA_FILE = '/Users/freud/Documents/data/PANC_TCGA/DATA/DrugSensitive_for_plot.csv'
#DATA_FILE = '/Users/freud/Documents/data/PANC_TCGA/DATA/DrugSensitive_log2_for_plot.csv'
#CALC_COL = ['TumorPurity']
#CALC_COL = ['StromalScore', 'ImmuneScore', 'ESTIMATEScore']
CALC_COL = list(pd.read_csv(DATA_FILE, index_col=[0]).columns[: -3])
#CALC_COL = list(pd.read_csv(DATA_FILE, index_col=[0]))

CLASS_FILE = '/Users/freud/Documents/data/PANC_TCGA/RESULT/LASSO_CLASS_RESULT/LASSO_CLASS_COMB.csv'
CLASS_COL = 'COMB'

dataframe = pd.read_csv(DATA_FILE)
classframe = pd.read_csv(CLASS_FILE)
classframe = classframe[['Unnamed: 0', CLASS_COL]]
classframe.columns = ['case_id', 'classname']

dataframe = pd.merge(dataframe, classframe, on='case_id', how='inner')
plateframe = Platedata(dataframe)

plt = DrawMutiBoxPlot(plateframe,
                      statcols=['hue', 'classname'],
                      valuecol='value',
                      statitem=CALC_COL)
plt.xticks(rotation=90)
plt.tight_layout()
plt.savefig(CFG.lassoresultpath + CLASS_COL + '_' + CALC_COL[0] + '.jpg')
plt.close()
