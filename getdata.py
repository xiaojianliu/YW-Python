# -*- coding: utf-8 -*-
"""
routines to get data 
"""
from matplotlib.dates import num2date
from pydap.client import open_url
import time
import numpy as np
import sys
import datetime as dt
import scipy

def getemolt_latlon(site):
    """
    get data from emolt_sensor 
    """
    urllatlon = 'http://gisweb.wh.whoi.edu:8080/dods/whoi/emolt_site?emolt_site.SITE,emolt_site.LAT_DDMM,emolt_site.LON_DDMM,emolt_site.ORIGINAL_NAME,emolt_site.BTM_DEPTH&emolt_site.SITE='
    dataset = open_url(urllatlon+'"'+site+'"')
    print dataset
    var = dataset['emolt_site']
    lat = list(var.LAT_DDMM)
    lon = list(var.LON_DDMM)
    original_name = list(var.ORIGINAL_NAME)
    bd=list(var.BTM_DEPTH)  
    return lat[0], lon[0], original_name,bd
    
def getemolt_temp(site):#,input_time=[dt.datetime(1880,1,1),dt.datetime(2020,1,1)],dep=[0,1000]):#you can input two values:start_time & end_time or one value:interval_days
  #gets eMOLT data using pydap, return datetime, temp, u, v, depth
  #realized in Nov 2012 that I need to fetch data based on time as well as site but I haven't implemented that yet
  url='http://gisweb.wh.whoi.edu:8080/dods/whoi/emolt_sensor?emolt_sensor.SITE,emolt_sensor.YRDAY0_LOCAL,emolt_sensor.TIME_LOCAL,emolt_sensor.TEMP,emolt_sensor.DEPTH_I,emolt_sensor.U,emolt_sensor.V&emolt_sensor.SITE='
  # get the emolt_sensor data
  dataset=open_url(url+'"'+site+'"')
  var=dataset['emolt_sensor']
  print 'hold on  ... extracting your eMOLT mooring data'
  temp=list(var.TEMP)
  #depth=list(var.DEPTH_I)
  time0=list(var.YRDAY0_LOCAL)
  year_month_day=list(var.TIME_LOCAL)
  print 'now generating a datetime'
  date_time=[]
  for i in range(len(time0)):
    date_time.append(num2date(time0[i]+1.0).replace(year=time.strptime(year_month_day[i],'%Y-%m-%d').tm_year).replace(month=time.strptime(year_month_day[i],'%Y-%m-%d').tm_mon).replace(day=time.strptime(year_month_day[i],'%Y-%m-%d').tm_mday))
  print 'sorting eMOLT time..'
  inds = np.argsort(date_time)
  date_time = np.take(date_time,inds)
  temp=np.take(temp,inds)
  return date_time,temp#,depth 
  
def getemolt_data(site):
    url='http://gisweb.wh.whoi.edu:8080/dods/whoi/emolt_sensor?emolt_sensor.SITE,emolt_sensor.YRDAY0_LOCAL,emolt_sensor.TIME_LOCAL,emolt_sensor.TEMP,emolt_sensor.DEPTH_I,emolt_sensor.U,emolt_sensor.V&emolt_sensor.SITE='
  # get the emolt_sensor data
    dataset=open_url(url+'"'+site+'"')
    var=dataset['emolt_sensor']
    print 'hold on  ... extracting your eMOLT mooring data'
    temp=list(var.TEMP)
    #depth=list(var.DEPTH_I)
    time0=list(var.YRDAY0_LOCAL)
    year_month_day = list(var.TIME_LOCAL)
    print 'now generating a datetime'
    return temp,time0,year_month_day  
    
def get_dataset(url):
    try:
        dataset = open_url(url)
    except:
        print 'Sorry, ' + url + 'is not available' 
        sys.exit(0)
    return dataset

def getdrift(id):
    """
    uses pydap to get remotely stored drifter data given an id number
    """
    url = 'http://gisweb.wh.whoi.edu:8080/dods/whoi/drift_data'
    
    dataset = get_dataset(url) 
     
    try:
        lat = list(dataset.drift_data[dataset.drift_data.ID == id].LAT_DD)
    except:
        print 'Sorry, ' + id + ' is not available'
        sys.exit(0)
        
    lon = list(dataset.drift_data[dataset.drift_data.ID == id].LON_DD)
    time0 = list(dataset.drift_data[dataset.drift_data.ID == id].TIME_GMT)
    yearday = list(dataset.drift_data[dataset.drift_data.ID == id].YRDAY0_GMT)
    dep = list(dataset.drift_data[dataset.drift_data.ID == id].DEPTH_I)
    datet = []
    # use time0('%Y-%m-%d) and yearday to calculate the datetime
    for i in range(len(yearday)):      
        datet.append(num2date(yearday[i]).replace(year=time.strptime(time0[i], '%Y-%m-%d').tm_year).replace(day=time.strptime(time0[i], '%Y-%m-%d').tm_mday))
    
    return lat, lon, datet, dep, time0, yearday  
        
    
def getobs_tempsalt(site,input_time,dep):
    """
    get data from url, return datetime, temperature, and start and end times
    input_time can either contain two values: start_time & end_time OR one value:interval_days
    """
    url = 'http://gisweb.wh.whoi.edu:8080/dods/whoi/emolt_sensor?emolt_sensor.SITE,emolt_sensor.YRDAY0_LOCAL,emolt_sensor.TIME_LOCAL,emolt_sensor.TEMP,emolt_sensor.DEPTH_I,emolt_sensor.SALT&emolt_sensor.SITE='
    dataset = get_dataset(url + '"' + site + '"')
    var = dataset['emolt_sensor']
    print 'extracting eMOLT data using PyDap... hold on'
    temp = list(var.TEMP)
    depth = list(var.DEPTH_I)
    time0 = list(var.YRDAY0_LOCAL)
    year_month_day = list(var.TIME_LOCAL)
    salt=list(var.SALT)
  
    print 'Generating a datetime ... hold on'
    datet = []
    for i in scipy.arange(len(time0)):
        #datet.append(num2date(time0[i]+1.0).replace(year=time.strptime(year_month_day[i], '%Y-%m-%d').tm_year).replace(day=time.strptime(year_month_day[i], '%Y-%m-%d').tm_mday))
        datet.append(num2date(time0[i]+1.0).replace(year=dt.datetime.strptime(year_month_day[i], '%Y-%m-%d').year).replace(month=dt.datetime.strptime(year_month_day[i],'%Y-%m-%d').month).replace(day=dt.datetime.strptime(year_month_day[i],'%Y-%m-%d').day).replace(tzinfo=None))
    #get the index of sorted date_time
    print 'Sorting mooring data by time'
    index = range(len(datet))
    index.sort(lambda x, y:cmp(datet[x], datet[y]))
    #reorder the list of date_time,u,v
    datet = [datet[i] for i in index]
    temp = [temp[i] for i in index]
    depth = [depth[i] for i in index]
    salt=[salt[i] for i in index]

    print 'Delimiting mooring data according to user-specified time'  
    part_t,part_time,part_salt = [], [],[]
    if len(input_time) == 2:
        start_time = input_time[0]
        end_time = input_time[1]
    if len(input_time) == 1:
        start_time = datet[0]
        end_time = start_time + input_time[0]
    if  len(input_time) == 0:
        start_time = datet[0]
        end_time=datet[-1]
    print datet[0], datet[-1]
    for i in range(0, len(temp)):
        if (start_time <= datet[i] <= end_time) & (dep[0]<=depth[i]<= dep[1]):
            part_t.append(temp[i])
            part_time.append(datet[i]) 
            part_salt.append(salt[i])
    temp=part_t
    datet=part_time
    salt=part_salt
    return datet,temp,salt