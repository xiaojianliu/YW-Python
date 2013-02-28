# -*- coding: utf-8 -*-
"""
Created on Thu Feb 21 13:24:56 2013

@author: jmanning
"""

import pandas as pd
from pandas.core.common import save
site=['AB01','BA01','BA02','BA03','BF01','BI01','BM01','BN01','CP01','DC01','DJ01']
df1=pd.read_csv(site[0]+'_wtmp_mc_mod_mc_obs.csv',index_col=0)
dfmean1=df1.mean()
df2=df1['mean'].fillna(0)
dfmean1['max']=max(df2)
dfmean1['min']=min(df2)
pdmean1=pd.DataFrame(dfmean1)
pdmean1=pdmean1.drop('median')
pdmean1.columns=[site[0]]

for k in range(len(site)):
    if k!=0:
        df=pd.read_csv(site[k]+'_wtmp_mc_mod_mc_obs.csv',index_col=0)
        dfmean=df.mean()
        df3=df['mean'].fillna(0)
        dfmean['max']=max(df3)
        dfmean['min']=min(df3)
        pdmean=pd.DataFrame(dfmean)
        pdmean=pdmean.drop('median')
        pdmean.columns=[site[k]]
        pdmean1=pdmean1.join(pdmean)
print pdmean1
pdmean1.to_csv('totalcaculate.csv',index=True)
htmlmean=pdmean1.to_html(header=True,index=True,float_format=lambda x: '%10.2f' % x)
save(htmlmean,'/net/nwebserver/epd/ocean/MainPage/lob/emoltvsnecofs_diff.html')