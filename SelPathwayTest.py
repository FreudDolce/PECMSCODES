#!/usr/bin/env python
# -*- coding: utf-8 -*-
# @Author : Ji Hongchen
# @Email : jhca.xyt@163.com
# @Last version : 2021-11-02 16:08
# @Filename : SelPathwayTest.py

import RandomSearchpathway
from RandomSearchpathway import GenerateRandomResult
import argparse
import cfg
import time

CFG = cfg.cfg()

parser = argparse.ArgumentParser()
parser.add_argument('-g', help='Gmt file name.')
args = parser.parse_args()

t = time.localtime()
now = str(t.tm_year) + '_' + str(t.tm_mon) + '_' + \
    str(t.tm_mday) + '_' + str(t.tm_hour) + '_' + str(t.tm_min)

GenerateRandomResult(args.g, args.g, now)
