"""
9/JAN/2013
Description:this is to read-in data and plot 'time-speed,direction'figure.

Import problem:import'getdata'
               replace'from getdata import gomoos'as'from getdata import *'
               
Response:no response since no 'pydap'

"""
import sys
# add the path of
sys.path.append("c:/oceanography")
import scipy
from getdata import *
from conversions import uv2sd
import matplotlib as mpl

interval_day=28
#per_panel_days=14
#site_str=[]
#site_list=[1,2,3,5,6,7,27, 28, 29, 30, 31, 32]
site_str=['A0125']
#for i in site_list:
  #  if i <10: i="0"+str(i)
   # site_str.append('"'+'NF'+str(i)+'"')
for site in site_str:
  if site=='A0125':
   # (time_obs,u_obs,v_obs,depth,sdate,edate)=getemolt(site,interval_day)
    #(lat,lon,original_name)=getemolt_latlon(site)
    #depth_i=-(np.mean(depth))
    #(mlat,mlon)=dm2dd(lat,lon)
    #print lat
 #  getgomoos
    (time_getgomoos_obs,u_getgomoos_obs,v_getgomoos_obs,mlat,mlon,sdate,edate,depth_i)=getgomoos(site,interval_day)
    depth_ii=-3
    speed,direction=[],[]
    for i in range(len(u_getgomoos_obs)):
       s=uv2sd(u_getgomoos_obs[i],v_getgomoos_obs[i])[0]
       d=uv2sd(u_getgomoos_obs[i],v_getgomoos_obs[i])[1]
       speed.append(s)
       direction.append(d)
   
    fig = plt.figure(1)
    ax1 = fig.add_subplot(111)
    ax1.plot(time_getgomoos_obs,speed,"-",label="speed")
    ax1.plot(time_getgomoos_obs,direction,"-",label="direction")
    locator = mpl.dates.AutoDateLocator()
    ax1.xaxis.set_major_locator(locator)
    monthsFmt = DateFormatter('%m/%d')
    ax1.xaxis.set_major_formatter(monthsFmt)
    ax1.axis(xmin=time_getgomoos_obs[0],xmax=time_getgomoos_obs[len(time_getgomoos_obs)-1])
    plt.legend()
    plt.show()
