# -*- coding: utf-8 -*-
"""
Created on Wed Feb 20 09:15:34 2013
documented at http://www.nefsc.noaa.gov/epd/ocean/MainPage/py/modvsobs/html/index.html 
@author: jmanning
"""


import matplotlib.pyplot as plt
from getdata import getemolt_latlon,getobs_tempsalt
from conversions import dm2dd,f2c
from utilities import my_x_axis_format
import pandas as pd
from numpy import float64
from datetime import datetime, timedelta
from models import getFVCOM_bottom_tempsalt_netcdf
import datetime as dt
import numpy as np
import netCDF4
def nearlonlat(lon,lat,lonp,latp):
    """
    i,min_dist=nearlonlat(lon,lat,lonp,latp)
    find the closest node in the array (lon,lat) to a point (lonp,latp)
    input:
        lon,lat - np.arrays of the grid nodes, spherical coordinates, degrees
        lonp,latp - point on a sphere
        output:
            i - index of the closest node
            min_dist - the distance to the closest node, degrees
            For coordinates on a plane use function nearxy
            
            Vitalii Sheremet, FATE Project
    """
    cp=np.cos(latp*np.pi/180.)
    # approximation for small distance
    dx=(lon-lonp)*cp
    dy=lat-latp
    dist2=dx*dx+dy*dy
    # dist1=np.abs(dx)+np.abs(dy)
    i=np.argmin(dist2)
    #min_dist=np.sqrt(dist2[i])
    return i#,min_dist 


def resamda(oritso):
        '''
        resample daily average data
        add some new columns like 'count mean,median,max,min,std,yy,mm,dd'
        '''
        oritsoda=pd.DataFrame(np.square(oritso[0].values),index=oritso.index)
        #oritso['square']=np.square(oritso[0].values)
        rms1=oritsoda[0].resample('D',how=['count','sum'])
        rms=np.sqrt(rms1['sum']/rms1['count'])
        resamda=float64(oritso[0]).resample('D',how=['count','mean','median','min','max','std'])
        resamda['rms']=rms.values
        resamda.ix[resamda['count']<minnumperday,['mean','median','min','max','std','rms']] = 'NaN'
        resamda['yy']=resamda.index.year
        resamda['mm']=resamda.index.month
        resamda['dd']=resamda.index.day
        output_fmt=['yy','mm','dd','count','mean','median','min','max','std','rms']
        resamdaf=resamda.reindex(columns=output_fmt)  
        return resamdaf

def resamma(oritso):
        '''
        resample month average data
        '''
        oritsoma=pd.DataFrame(np.square(oritso[0].values),index=oritso.index)
        rms1=oritsoma[0].resample('m',how=['count','sum'])
        rms=np.sqrt(rms1['sum']/rms1['count'])
        resamma=float64(oritso[0]).resample('m',how=['count','mean','median','min','max','std'],kind='period')
        resamma['rms']=rms.values
        resamma.ix[resamma['count']<25*numperday,['mean','median','min','max','std','rms']] = 'NaN'
        resamma['yy']=resamma.index.year
        resamma['mm']=resamma.index.month
        resamma['dd']=15
        output_fmt=['yy','mm','dd','count','mean','median','min','max','std','rms']
        resammaf=resamma.reindex(columns=output_fmt)# found I needed to generate a new dataframe to print in this order
        return resammaf

def resamdc(resamda):
        '''
        resample daily climatology
        '''
        newindex=[]
        for j in range(len(resamda)):    
                newindex.append(resamda['mean'].index[j].replace(year=2000)) # puts all observations in the same year
        repd=pd.DataFrame(resamda['mean'].values,index=newindex)
        
        repd['square']=np.square(repd.values)
        rms1=repd['square'].resample('D',how=['count','sum'])
        rms=np.sqrt(rms1['sum']/rms1['count'])
        resamdc=repd[0].resample('D',how=['count','mean','median','min','max','std'])    #add columns for custom date format
        resamdc['rms']=rms.values   
        resamdc['yy']=0
        resamdc['mm']=resamdc.index.month
        resamdc['dd']=resamdc.index.day
        output_fmt=['yy','mm','dd','count','mean','median','min','max','std','rms']
        resamdcf=resamdc.reindex(columns=output_fmt)# found I needed to generate a new dataframe to print in this order
        return resamdcf

def resammc(resamdc):
        '''
        resample month climatology
        '''
        resammc=resamdc['mean'].resample('m',how=['mean','median'],loffset=timedelta(days=-15))
        resammc['count']=0
        resammc['min']=0.
        resammc['max']=0.
        resammc['std']=0.
        resammc['rms']=0.
        recount=resamdc['count'].resample('m',how=['mean'],loffset=timedelta(days=-15)).values
        remi=resamdc['min'].resample('m',how=['mean'],loffset=timedelta(days=-15)).values
        rema=resamdc['max'].resample('m',how=['mean'],loffset=timedelta(days=-15)).values
        restd=resamdc['std'].resample('m',how=['mean'],loffset=timedelta(days=-15)).values
        rerms=resamdc['rms'].resample('m',how=['mean'],loffset=timedelta(days=-15)).values
        for kk in range(len(resammc)):
           resammc['count'].values[kk]=recount[kk]
           resammc['min'].values[kk]=remi[kk]
           resammc['max'].values[kk]=rema[kk]
           resammc['std'].values[kk]=restd[kk]
           resammc['rms'].values[kk]=rerms[kk]
        resammc['yy']=0
        resammc['mm']=resammc.index.month
        resammc['dd']=0
        output_fmt=['yy','mm','dd','count','mean','median','min','max','std','rms']
        resammcf=resammc.reindex(columns=output_fmt)# found I needed to generate a new dataframe to print in this order
        return resammcf

def diffdadc(diff):
        day=[]
        for i in range(len(diff)):
            day.append(str(diff.index[i].year)+'-'+str(diffda.index[i].month)+'-'+str(diffda.index[i].day))
        #daydadctime=pd.DataFrame(day,index=diff.index)
        #diff=diff.join(daydadctime)
        diff['date']=day
        difff=diff.reindex(columns=['date','mean','median','min','max','std','rms'])
        #diff.columns=['date','mean','median','min','max','std']
        return difff


site=['BF01']
siteprocess=[]
depthinfor=[]
minnumperday=18
numperday=24
intend_to='temp'##############notice intend_to can be 'temp'or'salinity'
surf_or_bott='bott'
for k in range(len(site)):
#################read-in obs data##################################
        print site[k]
        [lati,loni,on,bd]=getemolt_latlon(site[k]) # extracts lat/lon based on site code
        print bd
        [lati,loni]=dm2dd(lati,loni)#converts decimal-minutes to decimal degrees
        print lati,loni
        '''           
        [obs_dt,obs_temp]=getemolt_temp(site[k]) # extracts time series
        obs_dtindex=[]
        for kk in range(len(obs_temp)):
            obs_temp[kk]=f2c(obs_temp[kk]) # converts to Celcius
            obs_dtindex.append(datetime.strptime(str(obs_dt[kk])[:10],'%Y-%m-%d'))
        '''
        if surf_or_bott=='bott':
            dept=[bd[0]-0.25*bd[0],bd[0]+.25*bd[0]]
            #dept=[14,33]
        else:
            dept=[0,5]
        (obs_dt,obs_temp,obs_salt,distinct_dep)=getobs_tempsalt(site[k], input_time=[dt.datetime(1880,1,1),dt.datetime(2010,12,31)], dep=dept)
        depthinfor.append(site[k]+','+str(bd[0])+','+str(distinct_dep[0])+'\n')
        obs_dtindex=[]
        if intend_to=='temp':            
            for kk in range(len(obs_temp)):
                obs_temp[kk]=f2c(obs_temp[kk]) # converts to Celcius
                obs_dtindex.append(datetime.strptime(str(obs_dt[kk])[:10],'%Y-%m-%d'))
            obstso=pd.DataFrame(obs_temp,index=obs_dtindex)
        else:
            for kk in range(len(obs_salt)):
                obs_dtindex.append(datetime.strptime(str(obs_dt[kk])[:10],'%Y-%m-%d'))
            obstso=pd.DataFrame(obs_salt,index=obs_dtindex)   
        print 'obs Dataframe is ready'

##################generate resample DataFrame and putput file################################################
        reobsdaf=resamda(obstso)
        reobsmaf=resamma(obstso)
        reobsdcf=resamdc(reobsdaf)
        reobsmcf=resammc(reobsdcf)
        reobsdaf.to_csv(site[k]+surf_or_bott+intend_to+'_da_obs.csv',index=False,header=False,na_rep='NaN',float_format='%10.2f')
        reobsmaf.to_csv(site[k]+surf_or_bott+intend_to+'_ma_obs.csv',index=False,header=False,na_rep='NaN',float_format='%10.2f')
        reobsdcf.to_csv(site[k]+surf_or_bott+intend_to+'_dc_obs.csv',index=False,header=False,na_rep='NaN',float_format='%10.2f')
        reobsmcf.to_csv(site[k]+surf_or_bott+intend_to+'_mc_obs.csv',index=False,header=False,na_rep='NaN',float_format='%10.2f')
###################read-in mod data#################################

        starttime=obs_dt[0].replace(tzinfo=None)
        endtime=obs_dt[-1].replace(tzinfo=None)
        if surf_or_bott=='bott':
            layer=44
            vname=intend_to
            urlfvcom = 'http://www.smast.umassd.edu:8080/thredds/dodsC/fvcom/hindcasts/30yr_gom3'
            nc = netCDF4.Dataset(urlfvcom)
            nc.variables
            lat = nc.variables['lat'][:]
            lon = nc.variables['lon'][:]
            times = nc.variables['time']
            jd = netCDF4.num2date(times[:],times.units)
            var = nc.variables[vname]
            print 'Now find the coincide timestample'
            inode = nearlonlat(lon,lat,loni,lati)
            print inode
            modindex=netCDF4.date2index([starttime,endtime],times,select='nearest')
            modtso=pd.DataFrame(var[modindex[0]:modindex[1],layer,inode],index=jd[modindex[0]:modindex[1]])
#                     modtso=getFVCOM_bottom_tempsalt_netcdf(lati,loni,starttime,endtime,layer=44,vname=intend_to)
        else:
                     modtso=getFVCOM_bottom_tempsalt_netcdf(lati,loni,starttime,endtime,layer=0,vname=intend_to)  
        
        print 'now filter data to make daily obs and mod coincide.'
##############generate resample DataFrame and putput file#############
        remoddaf=resamda(modtso)
        remoddaf['mean'] += reobsdaf['mean']*0
        remoddaf['median'] += reobsdaf['median']*0
        remoddaf['min'] += reobsdaf['min']*0
        remoddaf['max'] += reobsdaf['max']*0
        remoddaf['std'] += reobsdaf['std']*0
        remoddaf['rms'] += reobsdaf['rms']*0
        print 'now filter data to make monthly obs and mod coincide.'
        remodmaf=resamma(modtso)
        remodmaf['mean'] += reobsmaf['mean']*0
        remodmaf['median'] += reobsmaf['median']*0   
        remodmaf['min'] += reobsmaf['min']*0
        remodmaf['max'] += reobsmaf['max']*0
        remodmaf['std'] += reobsmaf['std']*0
        remodmaf['rms'] += reobsmaf['rms']*0
#             if remodmaf.index.year[i] in reobsmaf.index.year and remodmaf.index.month[i] in reobsmaf.index.month:
        remoddcf=resamdc(remoddaf)
        remodmcf=resammc(remoddcf)
        remoddaf.to_csv(site[k]+surf_or_bott+intend_to+'_da_mod.csv',index=False,header=False,na_rep='NaN',float_format='%10.2f')
        remodmaf.to_csv(site[k]+surf_or_bott+intend_to+'_ma_mod.csv',index=False,header=False,na_rep='NaN',float_format='%10.2f')
        remoddcf.to_csv(site[k]+surf_or_bott+intend_to+'_dc_mod.csv',index=False,header=False,na_rep='NaN',float_format='%10.2f')
        remodmcf.to_csv(site[k]+surf_or_bott+intend_to+'_mc_mod.csv',index=False,header=False,na_rep='NaN',float_format='%10.2f')
        print 'now prepare plot'
          
            
##############plot da compare figure##################
        fig=plt.figure(figsize=(16,10))
        ax=fig.add_subplot(211)
        ax.plot_date(reobsdaf.index,reobsdaf['mean'],fmt='-')
        plt.grid()
        ax.plot_date(remoddaf.index,remoddaf['mean'],fmt='-',color='red')#bottom most value equals 44
        if intend_to=='temp':
                plt.ylabel('degree c')
        else:
                plt.ylabel('salinity')
        plt.title('eMOLT site '+site[k]+surf_or_bott+intend_to+' vs FVCOM ')
        plt.legend(['observed','modeled'],loc='best')
#        plt.show()
###############plot mc compare figure##################
        TimeDelta=reobsmcf.index[-1]-reobsmcf.index[0]          
        ax1 = fig.add_subplot(212)
        my_x_axis_format(ax1, TimeDelta)
        ax1.plot_date(reobsmcf.index,reobsmcf['mean'],fmt='-')
        plt.grid()
        ax1.plot_date(remodmcf.index,remodmcf['mean'],fmt='-',color='red')#bottom most value equals 44
        if intend_to=='temp':
                plt.ylabel('degree c')
        else:
                plt.ylabel('salinity')
        plt.title('eMOLT site '+site[k]+surf_or_bott+intend_to+' vs FVCOM ')
        plt.legend(['observed','modeled'],loc='best')
        plt.show()
        plt.savefig(site[k]+surf_or_bott+intend_to+'_mod_obs.png')
############calculate the different#######################
            #output_fmt=[0,'mean','median','min','max','std']
        diffmc=reobsmcf-remodmcf
        diffmc['month']=range(1,13)
            #month=pd.DataFrame(range(1,13),index=diffmc.index)
        diffmcf=diffmc.reindex(columns=['month','mean','median','min','max','std','rms'])
            #diffmc.columns=['month','mean','median','min','max','std']
        diffmcf.to_csv(site[k]+surf_or_bott+intend_to+'_mc_mod_mc_obs.csv',index=False,header=True,na_rep='NaN',float_format='%10.2f')

        diffma=reobsmaf-remodmaf
        date=[]
        for i in range(len(diffma)):
                date.append(str(diffma.index[i].year)+'-'+str(diffma.index[i].month))
        diffma['Year-Month']=date
            #datetimepd=pd.DataFrame(date,index=diffma.index)
            #diffma=diffma.join(datetimepd)
        diffmaf=diffma.reindex(columns=['Year-Month','mean','median','min','max','std','rms'])
            #diffma.columns=['Year-Month','mean','median','min','max','std']
        diffmaf.to_csv(site[k]+surf_or_bott+intend_to+'_ma_mod_ma_obs.csv',index=False,header=True,na_rep='NaN',float_format='%10.2f')

        diffda=reobsdaf-remoddaf
        diffdc=reobsdcf-remoddcf

        diffdaf=diffdadc(diffda)
        diffdaf.to_csv(site[k]+surf_or_bott+intend_to+'_da_mod_da_obs.csv',index=False,header=True,na_rep='NaN',float_format='%10.2f')
        diffdcf=diffdadc(diffdc)
        diffdcf.to_csv(site[k]+surf_or_bott+intend_to+'_dc_mod_dc_obs.csv',index=False,header=True,na_rep='NaN',float_format='%10.2f')
        siteprocess.append(site[k])
        
print 'processed sites include'+str(siteprocess)+'.'
import pickle
f = open('ProcessedSite.csv', 'w')
pickle.dump(siteprocess, f)
f.close()
d=open('Depthinformation.csv','w')
pickle.dump(depthinfor,d)
d.close()
        