#!/usr/bin/env python
# -*- coding: utf-8 -*-
# @Author : Ji Hongchen
# @Email : jhca.xyt@163.com
# @Last version : 2021-11-11 20:07
# @Filename : getclassevrexp.py

import pandas as pd
import numpy as np
import os
import cfg


def ShowPathClusterResult(c_n, p_n):
    print('>>>>>>>>  ', c_n, ' classes in consider.')
    gsvainfo = pd.read_csv('gsva_result.csv')
    class_file = 'ClusterResult/cluster_' + str(c_n) + '.csv'
    pclst_file = 'ClusterResult/p_cluster_' + str(p_n) + '.csv'
    classinfo = pd.read_csv(class_file)
    plistinfo = pd.read_csv(pclst_file)
    for i in range(p_n):
        plist = list(plistinfo['Unnamed: 0'][plistinfo['x'] == i + 1])
        print('=========================================')
        print(str(i + 1), ', pathways is as below:')
        print(plist)
        print('-----------------------------------------')
        for c in range(c_n):
            cexpmean = gsvainfo[gsvainfo['Unnamed: 0'].isin(plist)][
                    list(classinfo['Unnamed: 0'][classinfo['x'] == c + 1])].mean().sum()
            print('class: ', c + 1, ' total is ', len(classinfo[classinfo['x'] == c + 1]), ' mean is:')
            print(cexpmean)
