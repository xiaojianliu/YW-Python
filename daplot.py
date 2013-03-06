# -*- coding: utf-8 -*-
"""
Created on Tue Mar  5 14:30:42 2013

@author: jmanning
"""
import pandas as pd
import matplotlib.pyplot as plt
from datetime import datetime

def parse(datet):
        #print datet[0:10]
        dt=datetime.strptime(datet[0:10], '%Y %m %d')
        return dt
###########AB01 site###############
df=pd.read_csv('AB01_wtmp_da_mod.csv',sep=',',parse_dates={'datet':[0,1,2]},index_col='datet',date_parser=parse,names=['yy','mm','dd','count','mean','median','min','max','std'])
fig=plt.figure(figsize=(15,10))
ax=fig.add_subplot(311)
ax.plot(df.index,df['mean'].values)
df00=pd.read_csv('AB01_wtmp_da_obs.csv',sep=',',parse_dates={'datet':[0,1,2]},index_col='datet',date_parser=parse,names=['yy','mm','dd','count','mean','median','min','max','std'])
ax.plot(df00.index,df00['mean'].values)
ax.set_ylabel('Temperature (degC)')
ax.set_title('AB01 site temperature daily record')
ax.grid(True)
plt.legend(['modeled','observed'],loc='upper right', bbox_to_anchor=(1, 1.30),
          ncol=3, fancybox=True, shadow=True)
##########BN01 site#################
df1=pd.read_csv('BN01_wtmp_da_mod.csv',sep=',',parse_dates={'datet':[0,1,2]},index_col='datet',date_parser=parse,names=['yy','mm','dd','count','mean','median','min','max','std'])
ax1=fig.add_subplot(312)
ax1.plot(df1.index,df1['mean'].values)
df11=pd.read_csv('BN01_wtmp_da_obs.csv',sep=',',parse_dates={'datet':[0,1,2]},index_col='datet',date_parser=parse,names=['yy','mm','dd','count','mean','median','min','max','std'])
ax1.plot(df11.index,df11['mean'].values)
ax1.set_ylabel('Temperature (degC)')
ax1.set_title('BN01 site temperature daily record')
ax1.grid(True)
###########JS06 site################
df2=pd.read_csv('JS06_wtmp_da_mod.csv',sep=',',parse_dates={'datet':[0,1,2]},index_col='datet',date_parser=parse,names=['yy','mm','dd','count','mean','median','min','max','std'])
ax2=fig.add_subplot(313)
ax2.plot(df2.index,df2['mean'].values)
df22=pd.read_csv('JS06_wtmp_da_obs.csv',sep=',',parse_dates={'datet':[0,1,2]},index_col='datet',date_parser=parse,names=['yy','mm','dd','count','mean','median','min','max','std'])
ax2.plot(df22.index,df22['mean'].values)
ax2.set_ylabel('Temperature (degC)')
ax2.set_xlabel('Year')
ax2.set_title('JS06 site temperature daily record')
ax2.grid(True)
plt.show()
plt.savefig('three sites temperature plot.png')

