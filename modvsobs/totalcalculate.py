# -*- coding: utf-8 -*-
"""
Created on Wed Sep 25 10:27:37 2013

@author: jmanning
"""
import pandas as pd
from pandas.core.common import save
import numpy as np
from math import sqrt

site=['AB01','AG01','BA01','BA02','BC01','BD01','BF01','BI02','BI01','BM01','BM02','BN01','BS02','CJ01','CP01','DC01','DJ01','DK01','DMF1','ET01','GS01','JA01','JC01','JS06','JT04','KO01','MF02','MM01','MW01','NL01','PF01','PM02','PM03','PW01','RA01','RM02','RM04','SJ01','TA14','TA15','TS01']
column=["mean","min","max","rms","std"]
dffinal = pd.DataFrame(columns=column)
for k in range(len(site)):
        df=pd.read_csv("all_text_outputfile_new/"+site[k]+'botttemp_mc_mod_mc_obs.csv',index_col=0)
#        ts=pd.read_csv(site[k]+'bottsalinity_mc_obs.csv',index_col=0,names=['yy','mm','dd','count','mean','median','min','max','std'])
        dfmean=df.mean()
        df3=df['mean'].fillna(0)
        dfmean['max']=max(df3)
        dfmean['min']=min(df3)
        dfmean=dfmean.drop('rms')
        dfmean=dfmean.drop('std')
#        dfmean['std']=ts['std'].mean()
        std=df['mean'].std()
        rms=sqrt((np.sum(np.square(df3)))/(np.count_nonzero(df3)))
        row= pd.DataFrame([dict(rms=rms,std=std), ])
        row=row.T
        pdmean=pd.DataFrame(dfmean)
        pdmean=pdmean.append(row)
        pdmean=pdmean.drop('median')
        pdmean.columns=[site[k]]
        dffinal=dffinal.append(pdmean.T)

dffinal.to_csv('totalcaculate.csv',index=True)
htmlmean=dffinal.to_html(header=True,index=True,float_format=lambda x: '%10.2f' % x)
#save(htmlmean,'/net/nwebserver/epd/ocean/MainPage/lob/emoltvsnecofs_diff.html')
print dffinal