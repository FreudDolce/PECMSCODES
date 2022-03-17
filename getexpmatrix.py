#!/usr/bin/env python
# -*- coding: utf-8 -*-
# @Author : Ji Hongchen
# @Email : jhca.xyt@163.com
# @Last version : 2021-09-27 00:11
# @Filename : getexpmatrix.py

import numpy as np
import pandas as pd
import os
import json

PATH = '/Users/freud/Documents/data/PANC_TCGA/selectedGeneExpression/'
fl = os.listdir(PATH)
ENS_DICT = np.load(
    '/Users/freud/Documents/data/GENE_NAME/GENE_LIST/ensembl_gene.npy',
    allow_pickle=True
).item()
ENS_ENZ = np.load('/Users/freud/Documents/data/GENE_NAME/GENE_LIST/ens2enzid.npy',
                  allow_pickle=True).item()
indexlist = list(ENS_ENZ)

primrary_tumor = np.load(PATH + 'fileinfomation/file_list_primary_tumor.npy')
metastasis_tumor = np.load(PATH + 'fileinfomation/file_list_metastatic.npy')
normal_tissue = np.load(
    PATH + 'fileinfomation/file_list_solid_tissue_normal.npy')

kwargs = {'sep': '\t', 'comment': '#'}

jfile = '/Users/freud/Documents/data/PANC_TCGA/GeneExpression/files.2021-09-27.json'
TCGA_LIST = np.load('/Users/freud/Documents/data/PANC_TCGA/selectedGeneExpression/fileinfomation/TCGA_list.npy')
print(len(TCGA_LIST))


def GetPatientDict(jsonfile):
    filename2caseid = {}
    with open(jsonfile, 'r') as fp:
        jsondata = json.load(fp)
    for i in jsondata:
        filename2caseid[i['file_name']] = \
            i['cases'][0]['case_id']
    return filename2caseid


def SplitEns():
    for f in fl:
        if '.csv' in f:
            finfo = pd.read_csv(PATH + f)
            finfo['ENS'] = finfo['ENS'].str.split('.', expand=True)[0]
            finfo.to_csv(PATH + f.split('.')[0] + '.csv', index=False)
            print(f + ' finished')


def ReplaceENSToName(f):
    finfo = pd.read_csv(f)
    orderframe = finfo[finfo['ENS'].isin(indexlist)]
    orderframe = orderframe.reset_index(drop=True)
    for i in range(len(orderframe)):
        orderframe['ENS'][i] = ENS_ENZ[orderframe['ENS'][i]]
    return (orderframe)


if __name__ == '__main__':
    #ss = GetPatientDict(jfile)
    #np.save('file2caseid.npy', ss)
    # SplitEns()
    """
    for f in fl:
        if f in primrary_tumor:
            if (('.csv' in f) and (f in TCGA_LIST)):
                finfo = pd.read_csv(PATH + f)
                try:
                    totalinfo = pd.merge(totalinfo, finfo, how='inner', on='ENS')
                except NameError:
                    totalinfo = finfo
                print(totalinfo.shape)
    totalinfo.to_csv('PANC_MERGED_MRNA.csv', index=False)
    """
    repf = ReplaceENSToName('PANC_MERGED_MRNA.csv')
    repf = repf.drop_duplicates('ENS', keep='first')
    repf.to_csv('SELECTED_MRNA.csv', index=False)
