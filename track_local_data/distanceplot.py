# -*- coding: utf-8 -*-
"""
Created on Tue May  7 14:05:43 2013

@author: jmanning
"""

import matplotlib.pyplot as plt
import pandas as pd
from pylab import unique
import numpy as np

FList = np.genfromtxt('drift_data/Drift_list.csv',dtype=None,names=['FNs'],delimiter=',')
FNS=list(FList['FNs'])
#FNS=['ID_19886970.csv']
for kk in range(len(FNS)):
    dftotal=pd.DataFrame(columns=['distkm','modlat','modlon','reallat','reallon','time','timeday'])
    df=pd.read_csv(FNS[kk],sep=',',skiprows=1,index_col=0,names=['distkm','modlat','modlon','reallat','reallon','time']) 
    fig = plt.figure()
    timelist=df.time.values
    numlist = [int(x) for x in timelist]
    numlist=unique(numlist)
    if   len(numlist)%2<>0:
          numlistloop=len(numlist)/2
    else:
          numlistloop=(len(numlist)-1)/2
##################get the model data########################################################################################
    for i in range(numlistloop):
           t0=numlist[2*i]
           t1=numlist[2*i+2]
           dfsep = df[(df.time < t1) & (df.time >=t0)] 
           grouper = dfsep.groupby(dfsep.time < t0+1)
           before, after = grouper.get_group(True), grouper.get_group(False)
           before['timeday']=before.time.values%1
           after['timeday']=after.time.values%1+1
           dfsepfinal=before.append(after)
           plt.plot(dfsepfinal.timeday.values,dfsepfinal.distkm.values,'b')
           dfsepfinal['day']=['%10.1f' % (x) for x in dfsepfinal.timeday]
           print dfsepfinal
           dftotal=dftotal.append(dfsepfinal)
    plt.plot(dftotal.groupby(dftotal.day).distkm.mean().index,dftotal.groupby(dftotal.day).distkm.mean().values,'r',linewidth=3)
    plt.title(FNS[kk][0:-4]+'obs vs mod drift distance ')
    plt.savefig(FNS[kk][0:-4]+'distance.png')
    plt.xlabel('day')
    plt.ylabel('distance(km)')
    plt.show()
       