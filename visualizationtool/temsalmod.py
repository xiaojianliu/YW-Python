# -*- coding: utf-8 -*-
"""
Created on Thu May 23 09:10:22 2013

@author: jmanning
"""

from pylab import *
from matplotlib.collections import PolyCollection
import matplotlib.tri as Tri
from mpl_toolkits.basemap import Basemap
import datetime as dt
import netCDF4
import sys
import numpy as np

urlname=open("ctrl_temsalmod.csv", "r").readlines()[0][27:-1]
depth=int(open("ctrl_temsalmod.csv", "r").readlines()[1][22:-1])
TIME=open("ctrl_temsalmod.csv", "r").readlines()[2][31:-1]

if urlname=="30yr":
    
    stime=dt.datetime.strptime(TIME, "%Y-%m-%d %H:%M:%S")
    timesnum=stime.year-1981
    standardtime=dt.datetime.strptime(str(stime.year)+'-01-01 00:00:00', "%Y-%m-%d %H:%M:%S")
    timedeltaprocess=(stime-standardtime).days
    startrecord=26340+35112*(timesnum/4)+8772*(timesnum%4)+1+timedeltaprocess*24     
    url = 'http://www.smast.umassd.edu:8080/thredds/dodsC/fvcom/hindcasts/30yr_gom3?temp,lon,lat,lonc,latc,time,nv,h,siglay,salinity'
else:
    TIME=datetime.strptime(TIME, "%Y-%m-%d %H:%M:%S") 
    now=datetime.now()
    if TIME>now:
         diff=(TIME-now).days
    else:
         diff=(now-TIME).days
    if diff>3:
        print "please check your input start time,within 3 days both side form now on"
        sys.exit(0)
    numday=timedelta(days=numday)
    if TIME+numday>now+timedelta(days=3):
        print "please check your numday.access period is between [now-3days,now+3days]"
        sys.exit(0)
    timeperiod=(TIME+numday)-(now-timedelta(days=3))
    startrecord=(timeperiod.seconds)/60/60
    url="http://www.smast.umassd.edu:8080/thredds/dodsC/FVCOM/NECOFS/Forecasts/NECOFS_FVCOM_OCEAN_MASSBAY_FORECAST.nc?temp,lon,lat,lonc,latc,time,nv,h,siglay,salinity"
    
nc = netCDF4.Dataset(url)
lat = nc.variables['lat'][:]
lon = nc.variables['lon'][:]
latc = nc.variables['latc'][:]
lonc = nc.variables['lonc'][:]
temp=nc.variables['temp']
sali=nc.variables["salinity"]
siglay=nc.variables['siglay']
h = nc.variables['h'][:]
# read connectivity array
nv = nc.variables['nv'][:].T - 1
time_var = nc.variables['time']

# create a triangulation object, specifying the triangle connectivity array
tri = Tri.Triangulation(lon,lat, triangles=nv)
# plot depth using tricontourf
salinity=[]
temprature=[]
for i in range(len(lon)):
    depthtotal=siglay[:,i]*h[i]
    layer=np.argmin(abs(depthtotal+depth))
    print i,layer
    temprature.append(temp[startrecord,layer,i])
    sanility.append(sali[startrecord,layer,i])
temprature=np.array(temprature)
salinity=np.array(salinity)
fig=figure(figsize=(8,8))
ax=fig.add_subplot(211,aspect=1.0/cos(latc.mean() * pi / 180.0))
tricontourf(tri,temprature)
colorbar()
plt.title(urlname+' temp model') 
ax1=fig.add_subplot(212,aspect=1.0/cos(latc.mean() * pi / 180.0))
tricontourf(tri,temprature)
colorbar()
plt.title(urlname+' salinity model') 
plt.show()
