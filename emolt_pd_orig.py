# -*- coding: utf-8 -*-
"""
Created on Wed Dec 12 12:51:26 2012

@author: JiM
"""
import pandas as pd
import pylab
import matplotlib.pyplot as plt
from pandas import *
import numpy as np
from datetime import datetime,timedelta
from matplotlib.dates import num2date,date2num,DateFormatter
############read-in data######################
# HARDCODE
fn='mlc0102'
df = pd.read_csv('/net/home3/ocn/jmanning/py/yw/emolt/'+fn+'.txt', parse_dates={'datet':[0,1]}, skiprows=7,index_col='datet',keep_date_col=True,names=['Date','Time','Temp'])#Read in data as indexcolume is Datatime-index
plt.figure(figsize=(12,4))
plt.plot(df.index.to_pydatetime(),df['Temp'])#plot the figure by matplotlib
plt.title('Select approximate start and end points')
[start,end]=pylab.ginput(n=2)
temps=df.Temp
plt.close()

#####for start point zoom figure#########
startfront=temps[num2date(start[0]) - pd.offsets.Day(4):num2date(start[0]) + pd.offsets.Day(4)]#add the front 4 days' data and back 4 days's data
plt.plot(startfront.index,startfront[startfront.index])
plt.title("Click to select start time")
formatter = DateFormatter('%d-%b')
plt.gca().xaxis.set_major_formatter(formatter)
sfinal = pylab.ginput(1)
plt.close()
sfinaltime= num2date(sfinal[0][0])
print sfinaltime

#####for end point zoom figure###########
endback= temps[num2date(end[0]) - pd.offsets.Day(4):num2date(end[0]) + pd.offsets.Day(4)]
plt.plot(endback.index,endback[endback.index])
plt.title("Click to select end time")
plt.gca().xaxis.set_major_formatter(formatter)
efinal = pylab.ginput(1)
plt.close()
efinaltime = num2date(efinal[0][0])
print efinaltime

######for the final figure################
FF=df[sfinaltime:efinaltime]
criteria=3.0*FF.Temp.std() # standard deviations
a=0
for i in range(len(FF)-2):
    diff1=abs(FF.Temp[i+1]-FF.Temp[i]) 
    diff2=abs(FF.Temp[i+2]-FF.Temp[i+1])
    if diff1 > criteria and diff2 > criteria:
           print str(FF.index[i])+ ' is replaced by Nan'
           a+=1
           FF.Temp[i+1]=float('NaN')
print 'There are ' +str(a)+ ' points have replaced'
plt.figure(figsize=(12,4))
plt.plot(FF.index.to_pydatetime(),FF['Temp'])
plt.savefig(fn+'.png')
plt.savefig(fn+'.ps')
plt.show()
Ins='0'
Ps='0'
Sc='BJ01'
PFDATA={'sitecode':Series(Sc,index=FF.index),
        'instrument':Series(Ins,index=FF.index),
        'sernum':Series(fn[3:8],index=FF.index),#########should change according to different format############
        'probsetting':Series(Ps,index=FF.index),
        'Temp':Series(df['Temp'],index=FF.index),
        'Salinity':Series('99.999',index=FF.index)}
output_fmt=['sitecode','sernum','probsetting','instrument','Temp','Salinity']
PFDATADF=DataFrame(PFDATA)
PFDATADFNEW=PFDATADF.reindex(columns=output_fmt)        
PFDATADFNEW.to_csv(fn[0:-4]+'thisisdata.csv',float_format='%10.2f')
print "output file are ready"
