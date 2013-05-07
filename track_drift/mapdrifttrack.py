# -*- coding: utf-8 -*-
"""
Created on Mon Apr 15 09:36:40 2013

@author: jmanning
"""


import matplotlib.pyplot as plt
from mpl_toolkits.basemap import Basemap
import numpy as np
import pylab
from datetime import datetime
from pydap.client import open_url
#############get the index of lat and lon???
def nearlonlat(lon,lat,lonp,latp):
    cp=np.cos(latp*np.pi/180.)
    # approximation for small distance
    dx=(lon-lonp)*cp
    dy=lat-latp
    dist2=dx*dx+dy*dy
    #dist1=np.abs(dx)+np.abs(dy)
    i=np.argmin(dist2)
    #min_dist=np.sqrt(dist2[i])
    return i#,min_dist 
depth=-1
numdays=5
latsize=[40.5,44.5]
lonsize=[-71.,-67.0]
plt.figure(figsize=(7,6))
m = Basemap(projection='cyl',llcrnrlat=min(latsize)-0.01,urcrnrlat=max(latsize)+0.01,\
            llcrnrlon=min(lonsize)-0.01,urcrnrlon=max(lonsize)+0.01,resolution='h')#,fix_aspect=False)
m.drawparallels(np.arange(int(min(latsize)),int(max(latsize))+1,1),labels=[1,0,0,0])
m.drawmeridians(np.arange(int(min(lonsize)),int(max(lonsize))+1,1),labels=[0,0,0,1])
m.drawcoastlines()
m.fillcontinents(color='grey')
m.drawmapboundary()
plt.show()

###############################################################################
spoint = pylab.ginput(1)
la=spoint[0][1]
lo=spoint[0][0]
stime=datetime.strptime('2003-02-01 00:00:00', "%Y-%m-%d %H:%M:%S")
timesnum=stime.year-1981
standardtime=datetime.strptime(str(stime.year)+'-01-01 00:00:00', "%Y-%m-%d %H:%M:%S")
timedelta=(stime-standardtime).days
startrecord=26340+35112*(timesnum/4)+8772*(timesnum%4)+1+timedelta*24
endrecord=startrecord+24*numdays
latd,lond=[],[]

url='http://www.smast.umassd.edu:8080/thredds/dodsC/fvcom/hindcasts/30yr_gom3?'+'lon,lat,latc,lonc,Times['+str(startrecord)+':1:'+str(startrecord)+']'
#    url='http://www.smast.umassd.edu:8080/thredds/dodsC/fvcom/hindcasts/30yr_gom3?'+'h,lon,lat,latc,lonc,siglay,x,xc,y,yc,Times'+timeurl+',u'+timeurl+"[0:1:44][0:1:90414],"+"v"+timeurl+"[0:1:44][0:1:90414]"
dataset = open_url(url)
latc = np.array(dataset['latc'])
lonc = np.array(dataset['lonc'])
lat = np.array(dataset['lat'])
lon = np.array(dataset['lon'])
kf=nearlonlat(lonc,latc,lo,la) # nearest triangle center F - face
#kv=nearlonlat(lon,lat,lo,la)   # nearest triangle vertex V - vertex   
for i in range(startrecord,endrecord):
############read the particular time model from website#########    
    timeurl='['+str(i)+':1:'+str(i)+']'
    uvposition=str([0])+str([kf])
    url='http://www.smast.umassd.edu:8080/thredds/dodsC/fvcom/hindcasts/30yr_gom3?'+'Times'+timeurl+',u'+timeurl+uvposition+','+'v'+timeurl+uvposition   
    dataset = open_url(url)
    u=np.array(dataset['u'])
    v=np.array(dataset['v'])  
################get the point according the position###################
    print kf,u[0,0,0],v[0,0,0]
    par_u=u[0,0,0]
    par_v=v[0,0,0]   
    xdelta=par_u*60*60
    ydelta=par_v*60*60
    latdelta=ydelta/111111
    londelta=(xdelta/(111111*np.cos(la*np.pi/180)))
    la=la+latdelta
    lo=lo+londelta
    latd.append(la)
    lond.append(lo)
    kf=nearlonlat(lonc,latc,lo,la) # nearest triangle center F - face
#    kv=nearlonlat(lon,lat,lo,la)# nearest triangle vertex V - vertex  
m.plot(lond,latd,linewidth=1.5,color='r')
plt.title('model map track')
plt.show()
plt.savefig('driftrack.png')
