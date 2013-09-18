# -*- coding: utf-8 -*-
"""
Created on Fri Oct 12 15:27:09 2012

@author: vsheremet
"""
import numpy as np
#from pydap.client import open_url
import matplotlib.pyplot as plt
from SeaHorseLib import *
from datetime import *
#from scipy import interpolate
import sys
from SeaHorseTide import *
import shutil
import matplotlib.mlab as mlab
import matplotlib.cm as cm
from mpl_toolkits.basemap import Basemap
"""
from netCDF4 import Dataset

# read in etopo5 topography/bathymetry.
url = 'http://ferret.pmel.noaa.gov/thredds/dodsC/data/PMEL/etopo5.nc'
etopodata = Dataset(url)

topoin = etopodata.variables['ROSE'][:]
lons = etopodata.variables['ETOPO05_X'][:]
lats = etopodata.variables['ETOPO05_Y'][:]
# shift data so lons go from -180 to 180 instead of 20 to 380.
topoin,lons = shiftgrid(180.,topoin,lons,start=False)
"""



FN='binned_model.npz'
Z=np.load(FN) 

FNdrifter="binned_drifter.npz"
D=np.load(FNdrifter)
ub_meand=D['ub_mean']
vb_meand=D['vb_mean']

xb=Z['xb']
yb=Z['yb']
ub_mean=Z['ub_mean']
ub_median=Z['ub_median']
ub_std=Z['ub_std']
ub_num=Z['ub_num']
vb_mean=Z['vb_mean']
vb_median=Z['vb_median']
vb_std=Z['vb_std']
vb_num=Z['vb_num']
total_std=ub_std+vb_std
bigv=np.argwhere((ub_meand*ub_meand+vb_meand*vb_meand)>0.1250)#exclude big velocity
for j in range(len(bigv)):
    ub_meand[bigv[j][0],bigv[j][1]]=np.nan
    vb_meand[bigv[j][0],bigv[j][1]]=np.nan
ub_mean+=ub_meand*0
vb_mean+=vb_meand*0
less=np.argwhere(ub_num<=3)######exclude too less number
for i in range(len(less)):
    ub_mean[less[i][0],less[i][1]]=np.nan
    vb_mean[less[i][0],less[i][1]]=np.nan

Z.close()
#cmap = matplotlib.cm.jet
#cmap.set_bad('w',1.)
xxb,yyb = np.meshgrid(xb, yb)
cc=np.arange(-1.5,1.500001,0.05)
#cc=np.array([-1., -.75, -.5, -.25, -0.2, -.15, -.1, -0.05, 0., 0.05, .1, .15, .2, .25, .5, .75, 1.])
latsize=[40.5,45.5]
lonsize=[-72,-66]
plt.figure()
m = Basemap(projection='cyl',llcrnrlat=min(latsize),urcrnrlat=max(latsize),\
            llcrnrlon=min(lonsize),urcrnrlon=max(lonsize),resolution='h')#,fix_aspect=False)
#m.drawparallels(np.arange(-80,-55,20),labels=[1,1,0,0])
#m.drawmeridians(np.arange(34,48,5),labels=[0,0,0,1])
m.drawparallels(np.arange(int(min(latsize)),int(max(latsize))+1,2),labels=[1,0,0,0])
m.drawmeridians(np.arange(int(min(lonsize)),int(max(lonsize))+1,2),labels=[0,0,0,1])
m.drawcoastlines()
m.fillcontinents(color='grey')
m.drawmapboundary()
ub = np.ma.array(ub_mean, mask=np.isnan(ub_mean))
vb = np.ma.array(vb_mean, mask=np.isnan(vb_mean))
#Q=m.quiver(xxb,yyb,ub.T,vb.T,scale=2.5)
#plt.quiverkey(Q,0.7,0.09,0.5, '50cm/s', labelpos='W')
QQ=m.contourf(xxb,yyb,total_std.T)
b = plt.colorbar(QQ, orientation='vertical')
#plt.title('Model_derived mean')
plt.title('Model_derived_std mean')
#plt.plot(coast_lon,coast_lat,'b-')
#plt.axis([-72,-66,40.5,45.5])
plt.show()
