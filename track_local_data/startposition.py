# -*- coding: utf-8 -*-
"""
Created on Mon Aug  5 09:30:25 2013

@author: jmanning
"""

import pandas as pd
import numpy as np
from datetime import datetime
import matplotlib.pyplot as plt
from mpl_toolkits.basemap import Basemap
import matplotlib

FList = np.genfromtxt('drift_data/FList.csv',dtype=None,names=['FNs'],delimiter=',')
FNS=list(FList['FNs'])
xi = np.arange(-66,-72,-0.25)
yi = np.arange(40.5,45.5,0.25)
#firsttime=np.zeros((len(xi),len(yi)))
startlon=[]
startlat=[]
numbers=[]
lat=np.array([])
lon=np.array([])
for k in range(len(FNS)):
    print FNS[k]
    df=pd.read_csv('drift_data/'+FNS[k],index_col=0,skiprows=9,sep=',',names=['ID','TimeRD','TIME_GMT','YRDAY0_GMT','LON_DD','LAT_DD','TEMP','DEPTH'])
    print 'start time is '+df['TIME_GMT'].values[0],'end time is '+df['TIME_GMT'].values[-1]
    dtstart = datetime.strptime(df['TIME_GMT'].values[0], "%Y-%m-%d")
    dtend=datetime.strptime(df['TIME_GMT'].values[-1], "%Y-%m-%d")
    daynum=dtend-dtstart
    if daynum.days<7 or dtstart<datetime.strptime('2003-01-01', "%Y-%m-%d"):
        print dtstart
        print FNS[k]+' did not float more than 7 days or start time is before 2003'        
        k=+1
    else:
        lat=np.append(lat,df['LAT_DD'].values[0])
        lon=np.append(lon,df['LON_DD'].values[0])

plt.figure(figsize=(7,6))
m = Basemap(projection='cyl',llcrnrlat=40.5,urcrnrlat=45.5,\
            llcrnrlon=-72,urcrnrlon=-66,resolution='h')
m.drawparallels(np.arange(int(40.5),int(45.5),1),labels=[1,0,0,0])
m.drawmeridians(np.arange(int(-72),int(-66),1),labels=[0,0,0,1])
m.drawcoastlines()
m.fillcontinents(color='grey')
m.drawmapboundary()

x=lon
y=lat
ix=np.digitize(x,xi)
iy=np.digitize(y,yi)
for iix in range(0,len(xi)):
    for iiy in range(0,len(yi)):
            k,=np.where((ix==iix) & (iy==iiy))
            if len(k)<>0:
                print iix,iiy
                startlon.append(xi[iix])
                startlat.append(yi[iiy])
                numbers.append(len(k))
index = range(len(numbers))
index.sort(lambda x, y:cmp(numbers[x], numbers[y]))
numbers = [numbers[i] for i in index]
startlon = [startlon[i] for i in index]
startlat = [startlat[i] for i in index]                
                
                
                
                
legendcount_max=0
legendcount_min=0
legendcount_mid=0
bottomcurrentsitecount=0
for i in range(len(startlon)):
    if numbers[i]==min(numbers) and legendcount_min==0:
        Q=m.scatter(startlon[i],startlat[i],5*numbers[i],marker='o',color='red',label=str(numbers[i])+' drifter drops')
        legendcount_min+=1
        
    if numbers[i]==(max(numbers)+min(numbers))/2 and legendcount_mid==0:
        Q=m.scatter(startlon[i],startlat[i],5*numbers[i],marker='o',color='red',label=str(numbers[i])+' drifter drops')
        legendcount_mid+=1

    if numbers[i]==max(numbers) and legendcount_max==0:
        Q=m.scatter(startlon[i],startlat[i],5*numbers[i],marker='o',color='red',label=str(numbers[i])+' drifter drops')
        legendcount_max+=1
    else:
        Q=m.scatter(startlon[i],startlat[i],5*numbers[i],marker='o',color='red')
xxb,yyb = np.meshgrid(xi, yi)

dfsite=pd.read_csv('drift_data/deepcurrentsite1.csv',sep=',',names=['lat','lon'])
for ii in range(len(dfsite)):
    (a,b)=divmod(float(dfsite['lat'][ii]),100)   
    aa=int(a)
    bb=float(b)
    dfsite['lat'][ii]=aa+bb/60
    (c,d)=divmod(float(dfsite['lon'][ii]),100)
    cc=int(c)
    dd=float(d)
    dfsite['lon'][ii]=cc+(dd/60)
    if bottomcurrentsitecount==0:
        QQ=m.scatter(-dfsite['lon'][ii],dfsite['lat'][ii],20,marker='o',color='blue',label='bottom current site')
        bottomcurrentsitecount+=1
    else:
        QQ=m.scatter(-dfsite['lon'][ii],dfsite['lat'][ii],20,marker='o',color='blue')

dfsite1=pd.read_csv('drift_data/deepcurrentsite.csv',sep=',',names=['lat','lon'])
for iii in range(len(dfsite)):
    QQ=m.scatter(-dfsite['lon'][iii],dfsite['lat'][iii],20,marker='o',color='blue')

plt.title('Drifter start and current meter locations')
matplotlib.pyplot.legend(loc='lower right')
plt.show()
#plt.savefig(FNS[k][0:-4]+'.png')