# -*- coding: utf-8 -*-
"""
Created on Mon Jan 28 13:52:18 2013

@author: jmanning
"""
'''
step1:read-in data 
step2:groupby 'Year'and'Days'
step3:plot figure
'''
# emolt2_pd.py this is the Pandas code we use to generate plots for lobstermen
# Written by Yacheng Wang and Jim Manning
# This is a revision of the old "emolt2.py" which used scikits timeseries
from pydap.client import open_url
from pandas import DataFrame
from datetime import datetime as dt
from matplotlib.dates import num2date
import matplotlib.pyplot as plt
import sys
import matplotlib.dates as dates
import datetime
import pandas
from datetime import datetime
from dateutil import parser

#HARCODES####
site='BN01' # this is the 4-digit eMOLT site code you must know ahead of time
#############
def get_dataset(url):
    """
    Just checks to see if the OPeNDAP URL is working
    """
    try:
        dataset = open_url(url)
        print url+' is avaliable.'
    except:
        print 'Sorry, ' + url + 'is not available' 
        sys.exit(0)
    return dataset
###############
url = 'http://gisweb.wh.whoi.edu:8080/dods/whoi/emolt_sensor?emolt_sensor.SITE,emolt_sensor.YRDAY0_LOCAL,emolt_sensor.TIME_LOCAL,emolt_sensor.TEMP,emolt_sensor.DEPTH_I&emolt_sensor.SITE='
# get the emolt_sensor data
dataset = get_dataset(url + '"' + site + '"')
var = dataset['emolt_sensor']
print 'extracting eMOLT data using PyDap... hold on'
temp = list(var.TEMP)
depth = list(var.DEPTH_I)
time0 = list(var.YRDAY0_LOCAL)
year_month_day = list(var.TIME_LOCAL)
print 'Generating a datetime using the yearday records ... hold on'
datet = []
for i in range(len(time0)):
     datet.append(num2date(time0[i]+1.0).replace(year=dt.strptime(year_month_day[i], '%Y-%m-%d').year).replace(month=dt.strptime(year_month_day[i],'%Y-%m-%d').month).replace(day=dt.strptime(year_month_day[i],'%Y-%m-%d').day).replace(tzinfo=None))
tso=DataFrame(temp,index=datet)
#################add some keys for tso#########################
#tso['Month']=tso.index.month
tso['Year']=tso.index.year
tso['Day']=tso.index.dayofyear
#ax=tso.groupby(['Month','Year']).mean().unstack().plot(linewidth=3,legend=False)
tso1=tso.groupby(['Day','Year']).mean().unstack()
#################creat the datetime index#################################
date=[]
for i in range(len(tso1.index)-1):
#    date.append(datet[i].replace(year=2000))
    date.append(parser.parse(num2date(tso1.index[i]).replace(year=2000).isoformat(" ")))
date.append(parser.parse(num2date(tso1.index[len(tso1.index)-2]).replace(year=2000).isoformat(" ")))
'''
explain upside 
because tso1.index contain(1-366) so when we convert days to datetime format,366 will bacome 2000/jan/1,
so we delete the last index 366 and copy the last second record.
'''
######################################################################################################
ax = pandas.DataFrame(tso1.values,index=date).plot()
for i in range(len(ax.lines)):#plot in different ways
    if i<int(len(ax.lines)/2):
        ax.lines[i].set_linestyle('--')
        ax.lines[i].set_linewidth(2)
    elif i>=int(len(ax.lines)/2) and i<(len(ax.lines)-1):
        ax.lines[i].set_linestyle('-')
        ax.lines[i].set_linewidth(2)
    else:
        ax.lines[-1].set_linewidth(5)
        ax.lines[-1].set_color('black')
'''
below is to format the x-axis
'''        
ax.set_xlabel('Month')
#ax.title(site)
ax.set_title(site)
ax.xaxis.set_minor_locator(dates.MonthLocator(bymonth=None, bymonthday=1, interval=1, tz=None))
ax.xaxis.set_minor_formatter(dates.DateFormatter('%b'))
ax.xaxis.set_major_locator(dates.MonthLocator())
ax.xaxis.set_major_formatter(dates.DateFormatter(''))
patches,labels=ax.get_legend_handles_labels()
#ax.legend(set(tso['Year'].values),loc='center left', bbox_to_anchor=(.05, 0.5))
ax.legend(set(tso['Year'].values),loc='best')
plt.show()
plt.savefig('/net/nwebserver/epd/ocean/MainPage/lob/'+site+'.png')
plt.savefig('/net/nwebserver/epd/ocean/MainPage/lob/'+site+'.ps')