#!/usr/bin/env python
# -*- coding: utf-8 -*-
# @Author : Ji Hongchen
# @Email : jhca.xyt@163.com
# @Last version : 2021-10-03 20:24
# @Filename : ChangeClusterResult2Dict.py

import pandas as pd
import numpy as np
import os

workspace = '/Users/freud/Documents/data/PANC_TCGA/RESULT/'

fl = os.listdir(workspace)

for f in fl:
    fdict = {}
    if (('.csv' in f) and ('Cluster' in f)):
        cinfo = pd.read_csv(workspace + f)
        for i in cinfo.index:
            if cinfo['x'][i] in fdict:
                fdict[cinfo['x'][i]].append(cinfo['Unnamed: 0'][i])
            else:
                fdict[cinfo['x'][i]] = [cinfo['Unnamed: 0'][i]]
        np.save(f.split('.')[0] + '.npy', fdict)
