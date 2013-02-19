"""
8/JAN/2012
Description:read-in the wind data and plot 'time-u,v'figure then outputfile.

Import problem:NONE
"""
# Routine to plot NCEP reanalysis wind for a particular lat/lon and year
# Demonstrates the functions:
# "get_ncep_winds" in "c:/oceanography/models"
# "date2yd" in "c:/oceanography/conversions"
# "multiquiver" in "c:/oceanography/multiquiver"
# and writes to data to a file called "ncep_XXXX_wind.dat" where XXXX is the year
# huanxin January 26, 2012
# JiM modifications  Feb 10, 2012   
import sys
#sys.path.append("c:/oceanography")
import models
from models import get_ncep_winds
from conversions import *
from stickplot import multiquiver
import scipy

####### HARDCODES #######

yr=input('which year do you want?(2008):')
#yr=2008
syd=input('start yearday.(120):') #start yearday
#syd=120
eyd=input('end yearday.(250):') #end yearday
#eyd=250
lat=41.4
lon=-71
title='ncep_wind'
xlabel='time'
#########################
uw,vw,jdmat_w=get_ncep_winds(lat,lon,yr)
yeardays=date2yd(jdmat_w)
start_time=int(round(np.interp(syd,yeardays,range(len(yeardays)))))
end_time=int(round(np.interp(eyd,yeardays,range(len(yeardays)))))

u,v,time,y=[],[],[],[]      
for i in range(start_time,end_time):
    u.append(uw[i])
    v.append(vw[i])
    time.append(jdmat_w[i])
    y.append(0)
color="blue"
panels=1

plot=multiquiver(panels,time,u,v,color,title,xlabel)

#plot.savefig('jun_vs_jul_2009_wind.png')

#save_output_file('jun_vs_jul_2009_wind.dat', u, time, u, v, yeardays)
output_file='NCEP_'+str(yr)+'_wind.dat'
f=open(output_file,"w")
date_month,date_day,date_hour,date_second,date_minute=[],[],[],[],[]
for i in scipy.arange(len(time)):
  date_month.append(num2date(time[i]).month)
  date_day.append(num2date(time[i]).day)
  date_hour.append(num2date(time[i]).hour)
  date_second.append(num2date(time[i]).second)
  date_minute.append(num2date(time[i]).minute)
  f.write(str(date_month[i])+" "+str(date_day[i])+" " +
          str(date_hour[i])+" "+str(date_minute[i])+
          " "+str(u[i])+" "+str(v[i])+" "+"\n")
f.close()
