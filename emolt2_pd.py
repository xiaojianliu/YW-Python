# -*- coding: utf-8 -*-
"""
Pandas version of "emolt2"
Created on Thu Dec 13 10:51:36 2012

@author: jmanning & yacheng
"""

#imports
import pandas as pd
from pandas import *
from pylab import unique
import matplotlib.pyplot as plt
from pandas import Series
from matplotlib import mlab as ml
from matplotlib.dates import num2date,MonthLocator, DateFormatter
###################################################
# some hardcodes
hmday=7
numsamplestofilter =24
months    = MonthLocator()# every month
#monthsFmt = DateFormatter('%Y-%m-%d %H:%M:%S')
monthsFmt =DateFormatter('%b')
yrsFmt=DateFormatter('%Y')
output_fig_directory='/net/nwebserver/epd/ocean/MainPage/lob/'#the path of saving figure
####################################################
# load data from perl getts_emolt.plx site_id output
#T=ml.load('/net/home3/ocn/jmanning/atsee/this.dat') # file with year, yd,temp, depth
T=ml.load('this.dat')
#site=raw_input('Enter site code: ')
site='WD02'
surf_or_bot='bottom'#'surface (thin) and bottom (thick)'
dep=7.
dep2=3. # deliniates this process from geting other depths above dep2
print 'site='+site
####################################################

#generate timeseries from "T"
print 'Generating a datetime using the yearday records ... hold on'
datet = []
for i in range(len(T)):
    datet.append(num2date(T[i,1]+1.0).replace(year=int(T[i,0])))
tso=Series(T[:,2],index=datet)


# plotting all the years as if they were the SAME year
newindex=[]
for j in range(len(tso)):
    newindex.append(tso.index[j].replace(year=2000))

plt.figure()
with pd.plot_params.use('x_compat',True):

 for k in unique(tso.index.year):
    #id=ml.find(tso.index.year==k)
    t=Series(tso.values,index=newindex)
        
    print tso.index[0],tso.index[-1]
#t.plot(label=str(k))

plt.show()

#legend()    
