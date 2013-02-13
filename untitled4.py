# -*- coding: utf-8 -*-
"""
Created on Mon Jan 14 11:40:28 2013

@author: jmanning
"""

"""
Created on Mon Jan 14 10:27:06 2013

@author: jmanning
"""
import matplotlib as ml
import numpy as np
import matplotlib.pyplot as plt
import netCDF4
from scipy import spatial
import time

url = 'http://www.smast.umassd.edu:8080/thredds/dodsC/fvcom/hindcasts/30yr_gom3'

loni = -(70+((33+55/60.)/60.))
lati = 42+((31+23/60.)/60.)
nc = netCDF4.Dataset(url)
nc.variables

lat = nc.variables['lat'][:]
lon = nc.variables['lon'][:]
times = nc.variables['time']
jd = netCDF4.num2date(times[:],times.units)

# function to find index to nearest point
def near(array,value):
    idx=(np.abs(array-value)).argmin()
    return idx
    
def nearxy(y_array, x_array, y_point, x_point):
    distance = (y_array-y_point)**2 + (x_array-x_point)**2
    idy = np.where(distance==distance.min())
    return idy[0]

def do_kdtree(combined_x_y_arrays,points):
    mytree = spatial.cKDTree(combined_x_y_arrays)
    dist, indexes = mytree.query(points)
    return indexes

# find nearest point to desired location
inode = nearxy(lon,lat,loni,lati)

# Shoe-horn existing data for entry into KDTree routines
combined_x_y_arrays = np.dstack([lon.ravel(),lat.ravel()])[0]
points_list = [(loni,lati)]
inode2 = do_kdtree(combined_x_y_arrays,points_list)

# get all time records of variable [vname] at node [inode]
vname = 'temp'
var = nc.variables[vname]
time1=time.time()
#istop=20000
istop=len(jd)
chunk=10000
h=zeros(istop)
for i in range(0,istop,chunk):
    itime=range(i,min(i+chunk,istop))
    h[itime] = var[itime,0,inode]
    time2=time.time()
    print 'elapsed seconds = %.1f, %.1f percent done' % (time2-time1,(float(i)/istop*100.0)