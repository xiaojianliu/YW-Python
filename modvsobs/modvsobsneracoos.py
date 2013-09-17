# -*- coding: utf-8 -*-
"""
Created on Wed Feb 20 09:15:34 2013
documented at http://www.nefsc.noaa.gov/epd/ocean/MainPage/py/modvsobs/html/index.html 
@author: jmanning
"""






import pandas as pd
from numpy import float64
from datetime import datetime, timedelta
import numpy as np
import netCDF4
from utilities import nearlonlat

def resamda(oritso):
        '''
        resample daily average data
        add some new columns like 'count mean,median,max,min,std,yy,mm,dd'
        '''
        oritsoda=pd.DataFrame(np.square(oritso['salinity'].values),index=oritso.index)
        #oritso['square']=np.square(oritso[0].values)
        rms1=oritsoda[0].resample('D',how=['count','sum'])
        rms=np.sqrt(rms1['sum']/rms1['count'])
        resamda=float64(oritso['salinity']).resample('D',how=['count','mean','median','min','max','std'])
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
        oritsoma=pd.DataFrame(np.square(oritso['salinity'].values),index=oritso.index)
        rms1=oritsoma[0].resample('m',how=['count','sum'])
        rms=np.sqrt(rms1['sum']/rms1['count'])
        resamma=float64(oritso['salinity']).resample('m',how=['count','mean','median','min','max','std'],kind='period')
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

        
     
surf_or_bott="bott"
intend_to='salinity'        
minnumperday=18
numperday=24        
#site=['E02','E01',#'D02','B01','I01','F01','N01','M01']
site=['E01']
for k in range(len(site)):
#################read-in obs data##################################
        print site[k]
        def parse(datet):
            #print datet[0:10],datet[11:13],datet[14:16]
            dtime=datetime.strptime(datet[0:10],'%Y-%m-%d')
            delta=timedelta(hours=int(datet[11:13]),minutes=int(datet[14:16]))
            return dtime+delta
        obstso=pd.read_csv('sali_'+site[k]+'.csv',sep=',',skiprows=1,parse_dates={'datet':[0]},index_col='datet',date_parser=parse,names=['Date','depth','salinity','lat','lon'])
        print "NERACOOS observe data is ready! Then we try to resample observe data"
        reobsdaf=resamda(obstso)
        reobsmaf=resamma(obstso)
        reobsdcf=resamdc(reobsdaf)
        reobsmcf=resammc(reobsdcf)
        reobsdaf.to_csv(site[k]+surf_or_bott+intend_to+'_da_obs.csv',index=False,header=False,na_rep='NaN',float_format='%10.2f')
        reobsmaf.to_csv(site[k]+surf_or_bott+intend_to+'_ma_obs.csv',index=False,header=False,na_rep='NaN',float_format='%10.2f')
        reobsdcf.to_csv(site[k]+surf_or_bott+intend_to+'_dc_obs.csv',index=False,header=False,na_rep='NaN',float_format='%10.2f')
        reobsmcf.to_csv(site[k]+surf_or_bott+intend_to+'_mc_obs.csv',index=False,header=False,na_rep='NaN',float_format='%10.2f')
        print "now try to get data from internet"
        starttime=obstso.index[0].replace(tzinfo=None)
        endtime=obstso.index[-1].replace(tzinfo=None)
        lati=obstso.lat.values[0]
        loni=obstso.lon.values[0]
        depth=obstso.depth.values[0]
        
        
        urlfvcom = 'http://www.smast.umassd.edu:8080/thredds/dodsC/fvcom/hindcasts/30yr_gom3'
        nc = netCDF4.Dataset(urlfvcom)
        nc.variables
        lat = nc.variables['lat'][:]
        lon = nc.variables['lon'][:]
        times = nc.variables['time']
        siglay = nc.variables['siglay'][:]
        depthtotal=nc.variables['h'][:]
        jd = netCDF4.num2date(times[:],times.units)
        var = nc.variables[intend_to]
        print 'Now find the coincide timestample'
        inode = nearlonlat(lon,lat,loni,lati)
        deplay=(depthtotal[inode])*(siglay[:,inode])
        layer=np.argmin(abs(deplay+depth))
        modindex=netCDF4.date2index([starttime,endtime],times,select='nearest')
        modtso=pd.DataFrame(var[modindex[0]:modindex[1],layer,inode],index=jd[modindex[0]:modindex[1]]) 
        print "model fata frame was generated."         
        
        
        #modtso=getFVCOM_bottom_tempsalt_netcdf_neracoos(lati,loni,starttime,endtime,depth,intend_to)
        modtso.columns=['salinity']
        print 'now filter data to make daily obs and mod coincide.'
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
        remoddaf['rms'] += reobsdaf['rms']*0
        remoddcf=resamdc(remoddaf)
        remodmcf=resammc(remoddcf)
        remoddaf.to_csv(site[k]+surf_or_bott+intend_to+'_da_mod.csv',index=False,header=False,na_rep='NaN',float_format='%10.2f')
        remodmaf.to_csv(site[k]+surf_or_bott+intend_to+'_ma_mod.csv',index=False,header=False,na_rep='NaN',float_format='%10.2f')
        remoddcf.to_csv(site[k]+surf_or_bott+intend_to+'_dc_mod.csv',index=False,header=False,na_rep='NaN',float_format='%10.2f')
        remodmcf.to_csv(site[k]+surf_or_bott+intend_to+'_mc_mod.csv',index=False,header=False,na_rep='NaN',float_format='%10.2f')
        print 'now prepare plot'

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


            