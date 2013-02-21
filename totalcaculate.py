# -*- coding: utf-8 -*-
"""
Created on Thu Feb 21 13:24:56 2013

@author: jmanning
"""

import pandas as pd
from pandas.core.common import save
site=['AB01','BN01','AG01','BA01']
df1=pd.read_csv(site[0]+'_wtmp_mc_mod_mc_obs.csv',index_col=0)
dfmean1=df1.mean()
pdmean1=pd.DataFrame(dfmean1)
pdmean1.columns=[site[0]]
for k in range(len(site)):
    if k!=0:
        df=pd.read_csv(site[k]+'_wtmp_mc_mod_mc_obs.csv',index_col=0)
        dfmean=df.mean()
        pdmean=pd.DataFrame(dfmean)
        pdmean.columns=[site[k]]
        pdmean1=pdmean1.join(pdmean)
print pdmean1
pdmean1.to_csv('totalcaculate.csv',index=True)
htmlmean=pdmean1.to_html(header=True,index=True)
save(htmlmean,'totalcaculate.html')