# -*- coding: utf-8 -*-
"""
Created on Tue Feb 19 13:03:57 2013

@author: jmanning
"""

from pydap.client import open_url
from conversions import fth2m
from numpy import mean
import datetime as dt
from getdata import getemolt_temp,getemolt_latlon

site=['MM01','NARR','NF43','NL01','RB02','RM02','RM04','WD01','WHAQ']#,'DMF6']
for k in range(len(site)):
      [lat,lon,oriname]=getemolt_latlon(site,k)      
      [datet,temp,depth_i]=getemolt_temp(site,k,input_time=[dt.datetime(1880,1,1),dt.datetime(2020,1,1)], dep=[0,1000])
      depth=int(fth2m(mean(depth_i)))
      print '<tr><td style="vertical-align: top;"><a href='+site[k]+'_wtmp_da_'+str(depth)+'.csv>'+site[k]+'</a></td><td style="vertical-align: top;">Buzzards Bay</td><td style="vertical-align:top;"><a href='+site[k]+'_wtmp_dc_'+str(depth)+'.csv>'+str(lat)+'</td><td style="vertical-align: top;"><a href='+site[k]+'_wtmp_ma_'+str(int(depth))+'.csv>'+str(lon)+'</td><td style="vertical-align: top;"><a href='+site[k]+'_wtmp_mc_'+str(int(depth))+'.csv>'+str(int(depth))+'<br></td> </tr>' 
##################
#             <tr><td style="vertical-align: top;"><a href='DMF6_wtmp_da_25.csv'>DMF6</a></td><td style="vertical-align: top;">Buzzards Bay</td><td style="vertical-align: top;"><a href='DMF6_wtmp_dc_25.csv'>4126.79</td><td style="vertical-align: top;"><a href='DMF6_wtmp_ma_25.csv'>7059.35</td><td style="vertical-align: top;"><a href='DMF6_wtmp_mc_25.csv'>25<br></td> </tr>