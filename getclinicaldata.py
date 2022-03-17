#!/usr/bin/env python
# -*- coding: utf-8 -*-
# @Author : Ji Hongchen
# @Email : jhca.xyt@163.com
# @Last version : 2021-10-03 23:20
# @Filename : getclinicaldata.py

import pandas as pd
import numpy as np


clinc_info = ['case_id', 'case_submitter_id',
              'days_to_death', 'gender',
              'race', 'vital_status', 'year_of_death',
              'ajcc_pathologic_m', 'ajcc_pathologic_n', 'ajcc_pathologic_stage',
              'ajcc_pathologic_t',
              'days_to_last_follow_up', 'icd_10_code', 'morphology',
              'primary_diagnosis',
              'tissue_or_organ_of_origin']

kwargs = {'sep': '\t'}
clinical_info = pd.read_csv('/Users/freud/Documents/data/PANC_TCGA/selectedGeneExpression/fileinfomation/clinical.cart.2021-10-03/clinical.tsv', **kwargs)
selected_cinfo = clinical_info[clinc_info]
selected_cinfo = selected_cinfo.drop_duplicates(subset='case_id')
selected_cinfo.to_csv('clinical_info_for_R.csv', index=False)
