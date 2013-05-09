# -*- coding: utf-8 -*-
"""
Created on Thu May  9 12:33:11 2013

@author: jmanning
"""

import sys
sys.path.append("/home3/ocn/jmanning/py/jmanning/whython6/")
import basemap as bm
import matplotlib.pyplot as plt
import pylab
import numpy as np

plt.figure(4,figsize=(12,8))
bm.basemap_detail([40.0,43.0],[-71.,-65.0],True,True,0.5)
thismanager4 = pylab.get_current_fig_manager()
thismanager4.window.SetPosition((500, 500))
plt.title('basemap_detail')
plt.show()
#point=[]
#point.append(plt.ginput(n=20))
while True: 
    print "please click the points you want to save"
    n=raw_input("input the number of the points you want to save")
    point=plt.ginput(int(n))
    name=raw_input("input the bathymetry of the point you just pick up")
    np.save(name+'_depth',np.array(point))
    yesorno=raw_input("pick up again? y/n ")
    if yesorno == 'n':
           break
    if yesorno == 'y':
           continue
print "You have withdraw operation"
pointfile=["_70_depth.npy","_150_depth.npy"]
for kk in range(len(pointfile)):
     deppoint=np.load(pointfile[kk])
     lon=[]
     lat=[]
     for i in range(len(deppoint)):
          lon.append(deppoint[i][0])
          lat.append(deppoint[i][1])
     plt.plot(lon,lat,'ro')
     plt.show()
