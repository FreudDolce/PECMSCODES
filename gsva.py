#!/usr/bin/env python
# -*- coding: utf-8 -*-
# @Author : Ji Hongchen
# @Email : jhca.xyt@163.com
# @Last version : 2021-10-07 22:58
# @Filename : GSVA.py

import pandas as pd
import numpy as np
import os
import argparse
import cfg

CFG = cfg.cfg()

parser = argparse.ArgumentParser()
parser.add_argument('-g', help='Gmt file name.')
args = parser.parse_args()

expfile = CFG.EXPRESSION_FILE
gmtfile = CFG.GMT_PATH + args.g

os.system('Rscript ' + CFG.codepath + 'GSVA.R ' + expfile + ' ' + gmtfile + ' ' + CFG.resultpath)

gsva_result = pd.read_csv(CFG.resultpath + 'gsva_result.csv')
gsva_result.columns = gsva_result.columns
gsva_result.to_csv(CFG.resultpath + 'gsva_result.csv')

print('==================================================')
print('GSVA finished.')
print('==================================================')
