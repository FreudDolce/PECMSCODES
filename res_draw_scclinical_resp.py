#!/usr/bin/env python
# -*- coding: utf-8 -*-
# @Author : Ji Hongchen
# @Email : jhca.xyt@163.com
# @Last version : 2022-03-07 22:19
# @Filename : res_draw_scclinical_resp.py

import pandas as pd
import numpy as np
from res_clusterandclassdraw import DrawPileMap

scclinical_data = pd.read_csv('~/Documents/MANU/BI/PANC_ECM/DATA/Clinical_data_exp/Clinical_data.csv',
                              index_col=[0])

plt = DrawPileMap(scclinical_data,
                  class_col='PECMS_CLASS',
                  draw_col='BEST_RESPONSE',
                  draw_order=['PR', 'SD', 'PD'])

plt.legend(bbox_to_anchor=(0, 1.2))
plt.yticks(fontproperties='Arial', size=25)
plt.gcf().subplots_adjust(left=0.4, right=0.92, bottom=0.2)
plt.savefig(
    '/Users/freud/Documents/MANU/BI/PANC_ECM/Manuscript/Figure/SubFigure/scclinical_resp.tiff')
