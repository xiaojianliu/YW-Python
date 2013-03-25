# -*- coding: utf-8 -*-
"""
Created on Fri Mar 22 15:38:51 2013

@author: jmanning
"""

import pandas as pd
import matplotlib.pyplot as plt
from mpl_toolkits.basemap import Basemap
import numpy as np
######################################################   
dfcaculate=pd.read_csv('totalcaculate.csv',sep=',',skiprows=1,index_col=0,names=['site','mean','min','max','std'])
dfsite=pd.read_csv('site.csv',sep=',',index_col=0)
f = open('warmcold.csv', 'w')
for i in range(len(dfcaculate)):
    for j in range(len(dfsite)):
         if str(dfcaculate.index[i]) == str(dfsite.index[j]):
             f.write(dfcaculate.index[i]+','+str(dfsite[' emolt_site.LAT_DDMM'][j])+','+str(dfsite[' emolt_site.LON_DDMM'][j])+','+str(dfcaculate['mean'][i])+'\n')

f.close()
df=pd.read_csv('warmcold.csv',sep=',',skiprows=1,index_col=0,names=['site','lat','lon','mean'])
for i in range(len(df)):
    (a,b)=divmod(float(df['lat'][i]),100)   
    aa=int(a)
    bb=float(b)
    df['lat'][i]=aa+bb/60
    (c,d)=divmod(float(df['lon'][i]),100)
    cc=int(c)
    dd=float(d)
    df['lon'][i]=cc+(dd/60)
latsize=[40.0,45.0]
lonsize=[-72.,-67.0]
plt.figure(figsize=(7,6))
m = Basemap(projection='cyl',llcrnrlat=min(latsize)-0.01,urcrnrlat=max(latsize)+0.01,\
            llcrnrlon=min(lonsize)-0.01,urcrnrlon=max(lonsize)+0.01,resolution='h')#,fix_aspect=False)
m.drawparallels(np.arange(int(min(latsize)),int(max(latsize))+1,1),labels=[1,0,0,0])
m.drawmeridians(np.arange(int(min(lonsize)),int(max(lonsize))+1,1),labels=[0,0,0,1])
m.drawcoastlines()
m.fillcontinents(color='grey')
m.drawmapboundary()
x, y = m(-df['lon'].values,df['lat'].values)
for i in range(len(df)):
    if df['mean'][i]<=np.float64(0):
        print x[i],y[i] 
        m.scatter(x[i],y[i],30*df['mean'][i],marker='o',color='blue')
#    else:
#         m.scatter(x[i],y[i],30*df['mean'][i],marker='o',color='read')
    
plt.title('emolt temperature difference site')
plt.show()
plt.savefig('warmcold.png')