# -*- coding: utf-8 -*-
"""
Created on Mon Aug 12 10:00:08 2013

@author: jmanning
"""
import pylab
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.basemap import Basemap
import os
import math
import matplotlib as mpl
import numpy as np
import itertools
from pylab import *

xi = np.arange(-72,-66,0.25)
yi = np.arange(40.5,45.5,0.25)
xxb,yyb = np.meshgrid(xi, yi)

plt.figure(figsize=(7,6))
m = Basemap(projection='cyl',llcrnrlat=40.5,urcrnrlat=45.5,\
            llcrnrlon=-72,urcrnrlon=-66,resolution='h')
m.drawparallels(np.arange(int(40.5),int(46.5),0.25))#,labels=[1,0,0,0])
m.drawmeridians(np.arange(int(-72),int(-66),0.25))#,labels=[0,0,0,1])
m.drawparallels(np.arange(int(40.5),int(45.5),1),labels=[1,0,0,0])
m.drawmeridians(np.arange(int(-72),int(-66),1),labels=[0,0,0,1])
m.drawcoastlines()
m.fillcontinents(color='grey')
m.drawmapboundary()
for x in range(0,24):
    for y in range(0,20):
        if math.isnan(zb_mean[x,y])==False:
            plt.text(xxb[y,x]-0.25/2,yyb[y,x]+0.25/2,round(zb_mean[x,y],1),fontsize=8,fontweight='bold',
                    ha='center',va='center',color='b')
