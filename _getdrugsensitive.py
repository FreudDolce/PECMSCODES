#!/usr/bin/env python
# -*- coding: utf-8 -*-
# @Author : Ji Hongchen
# @Email : jhca.xyt@163.com
# @Last version : 2021-12-01 23:41
# @Filename : _getdrugsensitive.py

import os
import cfg
import pandas as pd
import numpy as np

CFG = cfg.cfg()

expdata = CFG.EXPRESSION_FILE

Drugs = ['Gemcitabine', 'Erlotinib', 'Docetaxel', 'Cisplatin', 'Olaparib',
         'Paclitaxel', '5-Fluorouracil', 'Thapsigargin', 'Phenformin']

for drug in Drugs:
    os.system('RScript ' + CFG.codepath + 'DRUGSENSITIVE.R ' +
              CFG.EXPRESSION_FILE + ' ' + drug + ' ' + CFG.datapath +
              'DrugSensitive/' + drug + '.csv')
