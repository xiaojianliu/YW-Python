# -*- coding: utf-8 -*-
"""
Created on Fri Feb 22 11:21:21 2013

@author: jmanning
"""


import pandas as pd
import matplotlib.pyplot as plt
from mpl_toolkits.basemap import Basemap
import numpy as np
######################################################   
threesite=['AB01','JS06','BN01']
df=pd.read_csv('ProcessedSite.csv',sep=',',skiprows=1,index_col=0,names=['site','lat','lon'])
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
m.scatter(x,y,20,marker='o',color='r')
for i in range(len(x)):
    if df.index[i] in threesite:
             print df.index[i]
             #plt.annotate(df.index[i],xy=(x[i],y[i]),xytext=(x[i]+0.51,y[i]-0.75),arrowprops=dict(facecolor='black',shrink=0.05));
             plt.annotate(df.index[i],xy=(x[i],y[i]),xytext=(x[i]+0.51,y[i]-0.75),arrowprops=dict(arrowstyle="-|>",connectionstyle="arc3,rad=-0.2",fc="w"));
if len(x)<=10:
     for i in range(len(x)):
         plt.text(x[i],y[i],df.index[i],fontsize=13,fontweight='normal',ha='right',va='baseline',color='b')
plt.title('emolt temperature site')
plt.show()
plt.savefig('EmoltSite.png')

   