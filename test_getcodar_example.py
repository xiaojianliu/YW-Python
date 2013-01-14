# -*- coding: utf-8 -*-
"""
Created on Fri Feb 10 11:45:04 2012

@author: jmanning
"""
"""
9/JAN/2013

Description:plot the 'lan-lon' figure 
    
Import problem:import getdata
    
Response:ImportError: No module named pydap.client

"""
# routine to track particles very crudely though the codar fields
# by extracting data from their respective OPeNDAP/THREDDS archives via the pyDAP method
#
# NOTE #1: There are previous versions of this code as it evolved as:
# "gettrack_JiM.py"

# NOTE #2: This code requires a few more python files that are included in the ocean.py including:
# -nearxy  to do 2-d interpolation
# -ll2uv  to caluclate u and v from lat, lon, and time
# Some non-essential calls can/should be commented out as follows:
#
# NOTE #4: There are several hardcodes for the user to change at the top of the code.
#
# Optional inputs include:
#
# drifterid  = distinct deployment code for drifter you are comparing with
#              where set to 0 if no drifter available     
#
# JiM       Oct 21 2011  

import sys
import os
import datetime
# add the path of ocean.py
import numpy
ops=os.defpath

if ops[2]<>'C':
    pydir='/net/home3/ocn/jmanning/'
else:
    pydir='c:/oceanography/'
sys.path.append(pydir)

from basemap import * 
from getdata import getcodar
import scipy

#### HARD CODE#################################
'''numdays=20 #number of days wanted to track
daystep=float(1)/24  #fraction of a day to step
x=1 #if you want give a box, and x is the number of the values in the box
'''
outputfile='x:/codar/'
'''# input the start of the time
startdate=date2num(datetime.datetime(2011,10,2))#,datetime.datetime(1981,10,15),datetime.datetime(1982,10,15),datetime.datetime(1983,10,15)]
print startdate
'''
###############################################  
#### input #########################################################
# choose_input=input("which do you want choose(you can input 1 or 2 or 3 or 4)? (1) one point, (2) list of points (3) give two points that define a box (4) give a box and depth")
la=40.0
lo=-73.0
datetime_wanted=datetime.datetime(2007,9,20)
ss=10 #number of vectors to skip or "subsample" 
                   
#for k in range(5):
#   datetime_wanted=datetime_wanted+datetime.timedelta(1) 
  # for different url 
# url_choose=input("which model? 0. CODAR 1. 30yr, 2. MASSBAY, 3. GOM3, 4. MASSBAY VS GOM3:")
url_choose=0
numhrs=2
if url_choose==0:
    url_name="CODAR"
    url="http://tashtego.marine.rutgers.edu:8080/thredds/dodsC/cool/codar/totals/macoora8km"
    #url="http://tashtego.marine.rutgers.edu:8080/thredds/dodsC/cool/codar/totals/macoora6km"
    #url="http://tashtego.marine.rutgers.edu:8080/thredds/dodsC/cool/codar/totals/sw06"
for k in range(numhrs):
  #datetime_wanted+=datetime.timedelta(0,3600)
     
  lat_vel,lon_vel,u,v=getcodar(url,datetime_wanted)
  datetime_wanted+=datetime.timedelta(1,21600)
  #print min(lat_vel),max(lat_vel)
  #print min(lon_vel),max(lon_vel)
  print type(lat_vel),numpy.size(lat_vel)
  #print len(lon_vel)
  print type(u),numpy.size(u)
  #print len(v)
  #print "past getcodar call with first time =+str(u[0][0])"
  
  #fig=plt.figure(1)
  fig=plt.subplot(3,2,k+1)
  idg=list(np.where(np.array(u)<>-999.0/100.))
  
  #print len(np.array(lon_vel)[idg[0]])
  #print len(np.array(lat_vel)[idg[0]])
  #print len(np.array(u)[idg[0]])
  #print len(np.array(v)[idg[0]])
  #lonv=list(np.array(lon_vel)[idg[0]])
  if len(idg)<>0:
    #print lon_vel[idg]
    
    q=plt.quiver(list(np.array(lon_vel)[idg]),list(np.array(lat_vel)[idg]),list(np.array(u)[idg]),list(np.array(v)[idg]),angles='xy',scale=10,color='r')
    plt.plot(lo,la,'k.',markersize=10)
    plt.title('') 
    plt.title('date & time: ' +str(datetime_wanted))
   
    bathy=False  

    #basemap_usgs([41,42],[-72,-76],bathy)
    #basemap_usgs([36,42],[-73,-76],bathy)
    #plt.plot()
  else:
    print "Sorry. No good data at this time"
    #lat_k,lon_k,time=gettrack_codar(jdmat_m,lon_vel,lat_vel,u,v,startdate,numdays,daystep,la,lo)
    #print "past gettrack_codar call"
    #plt.plot(lon_k,lat_k,"-",linewidth=1,marker=".",markerfacecolor='r')  
plt.show()
'''
    text1=plt.annotate(str(num2date(time[0]).month)+"-"+str(num2date(time[0]).day),
                                xy=(lon_k[0], lat_k[0]),  xycoords='data',
                                xytext=(-5, -8), textcoords='offset points',
                                arrowprops=dict(arrowstyle="->"))
    for i in range(1,len(time)):
         if num2date(time[i-1]).day<>num2date(time[i]).day: # annotates at midnight
            text1=plt.annotate(str(num2date(time[i]).month)+"-"+str(num2date(time[i]).day),          
                                        xy=(lon_k[i], lat_k[i]),  xycoords='data',
                                        xytext=(-8, -10), textcoords='offset points',
                                        arrowprops=dict(arrowstyle="->"))
          
        
    #set the numbers formatter
    majorFormatter1 = FormatStrFormatter('%.2f')
    majorFormatter2 = FormatStrFormatter('%.1f')
    #plt.yaxis.set_major_formatter(majorFormatter2)
    #plt.xaxis.set_major_formatter(majorFormatter1)
    plt.title(str(num2date(time[i-1]).year)+ " year " +url_name)
    #plt.xaxis.set_label_coords(0.5, -0.1)#set the position of the xlabel
    #plt.yaxis.set_label_coords(-0.1, 0.5)
'''    
#plt.savefig(outputfile+str(k)+'.png')
    #datetime_wanted+=datetime.timedelta(0,3600)
  
  
