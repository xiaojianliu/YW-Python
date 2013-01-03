# -*- coding: utf-8 -*-
"""
Created on Thu Dec  6 14:22:47 2012

@author: jmanning
"""
import pandas
import pylab
import matplotlib.pyplot as plt
from pandas import *
import numpy as np
from datetime import datetime,timedelta
from matplotlib.dates import num2date,date2num 
############read-in data######################
def parse(datet):
    #print datet[1:11],datet[14:19]
    dt=datetime.strptime(datet[1:11],'%Y-%m-%d')
    delta=timedelta(hours=int(datet[14:16]),minutes=int(datet[17:19]))
    return dt+delta
df=read_csv('c:/py/yacheng/mmf0401.txt',sep='\s+',skiprows=7,parse_dates={'datet':[0,1]},index_col='datet',date_parser=parse,names=['Date','Time','Temp'])
df.plot()
print "click on start and stop times to save date"
[start,end]=pylab.ginput(n=2)
print start,end

#####for start point zoom figure#########
startfront=start[0]-1
startback=start[0]+1
sfforplot=(num2date(startfront)).replace(minute=0,second=0,microsecond=0).isoformat(" ")
sbforplot=(num2date(startback)).replace(minute=0,second=0,microsecond=0).isoformat(" ")
df[sfforplot:sbforplot].plot()
sfinal=pylab.ginput(n=1)
print sfinal
sfinaltime=(num2date(sfinal[0][0])).replace(tzinfo=None)
sfinalforplot=sfinaltime.replace(minute=0,second=0,microsecond=0).isoformat(" ")
print sfinalforplot

#####for end point zoom figure###########
endfront=end[0]-1
endback=end[0]+1
efforplot=(num2date(endfront)).replace(minute=0,second=0,microsecond=0).isoformat(" ")
ebforplot=(num2date(endback)).replace(minute=0,second=0,microsecond=0).isoformat(" ")
df[efforplot:ebforplot].plot()
efinal=pylab.ginput(n=1)
print efinal
efinaltime=(num2date(efinal[0][0])).replace(tzinfo=None)
efinalforplot=efinaltime.replace(minute=0,second=0,microsecond=0).isoformat(" ")
print efinalforplot
plt.show()

