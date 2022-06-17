# -*- coding: utf-8 -*-
"""
Created on Tue Apr 26 09:21:49 2022
@author: Jean-Sébastien Gosselin
"""

from pyhelp.output import HelpOutput
import pandas as pd
import numpy as np

output_help = HelpOutput("D:/Projets/pyhelp/example/help_example.out")

rechg_monthly = output_help.data['rechg']

rechg_yearly = pd.DataFrame(
    data=np.sum(rechg_monthly, axis=2),
    index=output_help.data['cid'],
    columns=output_help.data['years']
    )
rechg_yearly.index.name = 'cid'

rechg_yearly.insert(0, 'lon_dd', output_help.data['lon_dd'])
rechg_yearly.insert(0, 'lat_dd', output_help.data['lat_dd'])

rechg_yearly.to_csv("D:/Projets/pyhelp/example/rechg_yearly.csv")
