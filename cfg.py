#!/usr/bin/env python
# -*- coding: utf-8 -*-
# @Author : Ji Hongchen
# @Email : jhca.xyt@163.com
# @Last version : 2021-10-10 20:12
# @Filename : cfg.py

class cfg():
    def __init__(self):
        self.worksapce = '/Users/freud/Documents/MANU/BI/PANC_ECM/'
        self.codepath = self.worksapce + 'code/'
        self.datapath = self.worksapce + 'DATA/'
        self.resultpath = self.worksapce + 'RESULT/'
        self.lassoresultpath = self.worksapce + 'RESULT/LASSO_CLASS_RESULT/'
        self.EXPRESSION_FILE = self.datapath + 'SELECTED_MRNA_SYMBOL.csv'
        self.CPT_EXP_FILE = self.datapath + 'CPTAC_MRNA_SYMBOL.csv'
        self.ORG_EXP_FILE = self.datapath + 'ORG_MRNA_SYMBOL.csv'
        self.SEARCHED_PATHWAY = self.datapath + 'Searched_pathway/'
        self.SELECTED_PATHWAY_PATH = self.datapath + 'SelectedPathway/'
        self.GMT_PATH = self.datapath + 'SelectedGmt/'
        self.ORI_GMT = self.datapath + \
            'gmt/msigdb_v7.4/msigdb_v7.4_GMTs/msigdb.v7.4.entrez.gmt'
        self.GENE_ID_TRANS_PATH = self.datapath + 'GENEIDTRANSFER/'
        self.CIBERSORT_RESULT = self.datapath + 'CIBERSORT_RESULT.csv'
        self.ESTIMATE_RESULT = self.datapath + 'PADC_estimate_score.csv'
        self.CLUSTER_RESULT_PATH = self.resultpath + 'ClusterResult/'
        self.MERGED_PATH = self.datapath + 'Merged_pathway/'
        self.clinicaldata= self.datapath + 'clinical_TCGA_for_R.csv'

        self.MIN_CLASS = 2
        self.MAX_CLASS = 4
        self.PATH_SELECT_PER_ITER = 12
        self.COL_THERS = 0.4876143877761896
        self.HA_THERS = 0.46174615286330845
        self.colors = ['#5372AB', '#D1875B', '#6AA56E', '#B65555', '#79AF97', '#6A6599', '#B0796B']
