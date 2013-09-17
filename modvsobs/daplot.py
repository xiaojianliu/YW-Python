# -*- coding: utf-8 -*-
"""
Created on Thu Mar  7 10:07:34 2013

@author: jmanning
"""
import pandas as pd
import matplotlib.pyplot as plt
from datetime import datetime

def parse(datet):
        #print datet[0:10]
        dt=datetime.strptime(datet[0:10], '%Y %m %d')
        return dt
###########read-in data as DataFrame###############
df=pd.read_csv("all_text_outputfile/E01bottsalinity_da_obs.csv",sep=',',parse_dates={'datet':[0,1,2]},index_col='datet',date_parser=parse,names=['yy','mm','dd','count','mean','median','min','max','std','rms'])
df00=pd.read_csv('all_text_outputfile/E01bottsalinity_da_mod.csv',sep=',',parse_dates={'datet':[0,1,2]},index_col='datet',date_parser=parse,names=['yy','mm','dd','count','mean','median','min','max','std','rms'])
#grouperdf00 = df00.groupby(df00.index < pd.Timestamp('2003-01-1'))
#beforedf00, afterdf00 = grouperdf00.get_group(True), grouperdf00.get_group(False)

#df1=pd.read_csv('BN01botttemp_da_mod.csv',sep=',',parse_dates={'datet':[0,1,2]},index_col='datet',date_parser=parse,names=['yy','mm','dd','count','mean','median','min','max','std'])
#df11=pd.read_csv('BN01botttemp_da_obs.csv',sep=',',parse_dates={'datet':[0,1,2]},index_col='datet',date_parser=parse,names=['yy','mm','dd','count','mean','median','min','max','std'])

#df2=pd.read_csv('all_text_outputfile/BF01botttemp_da_mod_da_obs.csv',sep=',',skiprows=1,index_col='date',parse_dates=True,names=['date','mean','median','min','max','std','rms'])
#df22=pd.read_csv('all_text_outputfile/I01botttemp_da_mod_da_obs.csv',sep=',',skiprows=1,index_col='date',parse_dates=True,names=['date','mean','median','min','max','std','rms'])
#grouperdf22 = df22.groupby(df22.index < pd.Timestamp('2003-01-1'))
#beforedf22, afterdf22 = grouperdf22.get_group(True), grouperdf22.get_group(False)

################plot the figure respectively########################

fig=plt.figure(figsize=(15,10))#figsize=(15,10))
ax=fig.add_subplot(111)
ax.plot(df.index,df['mean'].values,color='red')
ax.plot(df00.index,df00['mean'].values)
ax.set_ylabel('Salinity',fontsize=20)
ax.set_title('Bottom salinity at E01',fontsize=18)
ax.grid(True)
ax.lines[0].set_linewidth(3)
ax.lines[1].set_linewidth(2)
ax.tick_params(axis='both', which='major', labelsize=20)
plt.legend(['observed','modeled'],loc='upper right',# bbox_to_anchor=(1.01, 1.20),
          ncol=3, fancybox=True, shadow=True,prop={'size':14})
##########BN01 site#################
'''
ax1=fig.add_subplot(312,sharex=ax)
ax1.plot(df1.index,df1['mean'].values)
ax1.plot(df11.index,df11['mean'].values,color='red')
ax1.set_ylabel('Temperature (degC)')
ax1.set_title('Bottom Temperature at BN01')
ax1.grid(True)
for i in range(len(ax1.lines)):#plot in different ways
     ax1.lines[i].set_linewidth(2)

###########JS06 site################

ax2=fig.add_subplot(212,sharex=ax)
ax2.plot(df2.index,df2['mean'].values)
ax2.plot(afterdf22.index,afterdf22['mean'].values,color='red')
ax2.set_ylabel('temperature(degreeC)',fontsize=20)
ax2.set_xlabel('Year',fontsize=20)
ax2.set_title('I01(50m)vsBF01(60m) observed-modeled bottom temperature',fontsize=18)
ax2.grid(True)
#for i in range(len(ax2.lines)):#plot in different ways
#     ax2.lines[i].set_linewidth(2)
ax2.lines[0].set_linewidth(2)
ax2.lines[1].set_linewidth(2)
ax2.tick_params(axis='both', which='major', labelsize=20)
#plt.setp(ax1.get_xticklabels(),visible=False)
plt.setp(ax.get_xticklabels(),visible=False)
plt.tick_params(axis='both', which='major', labelsize=20)
'''
plt.show()
plt.savefig('E01bottsalinity.png')
