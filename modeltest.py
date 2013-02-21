# -*- coding: utf-8 -*-
"""
Created on Wed Feb 20 09:15:34 2013

@author: jmanning
"""
import matplotlib.pyplot as plt
import netCDF4
from getdata import getemolt_latlon,getemolt_temp
from conversions import dm2dd,f2c
from utilities import nearxy,my_x_axis_format
import pandas as pd
from numpy import float64
from datetime import datetime, timedelta





site=['AG01','BA01']
minnumperday=18
numperday=24
for k in range(len(site)):
        
        
        [lati,loni,on]=getemolt_latlon(site,k) # extracts lat/lon based on site code
        [lati,loni]=dm2dd(lati,loni) #converts decimal-minutes to decimal degrees
        [obs_dt,obs_temp]=getemolt_temp(site,k) # extracts time series
        for kk in range(len(obs_temp)):
            obs_temp[kk]=f2c(obs_temp[kk]) # converts to Celcius
        obs_dtindex=[]
        for i in range(len(obs_dt)):
            obs_dtindex.append(datetime.strptime(str(obs_dt[i])[:10],'%Y-%m-%d'))
        obstso=pd.DataFrame(obs_temp,index=obs_dtindex)
#########creat obs da##############
        reobsda=float64(obstso[0]).resample('D',how=['count','mean','median','min','max','std'])
        reobsda.ix[reobsda['count']<minnumperday,['mean','median','min','max','std']] = 'NaN'
        reobsda['yy']=reobsda.index.year
        reobsda['mm']=reobsda.index.month
        reobsda['dd']=reobsda.index.day
        output_fmt=['yy','mm','dd','count','mean','median','min','max','std']
        reobsdaf=reobsda.reindex(columns=output_fmt)  
        reobsdaf.to_csv(site[k]+'_wtmp_da_obs.csv',index=False,header=False,na_rep='NaN',float_format='%10.2f')
#########creat obs ma##############
        reobsma=float64(obstso[0]).resample('m',how=['count','mean','median','min','max','std'],kind='period')
        reobsma.ix[reobsma['count']<25*numperday,['mean','median','min','max','std']] = 'NaN'
        reobsma['yy']=reobsma.index.year
        reobsma['mm']=reobsma.index.month
        reobsma['dd']=15
        output_fmt=['yy','mm','dd','count','mean','median','min','max','std']
        reobsmaf=reobsma.reindex(columns=output_fmt)# found I needed to generate a new dataframe to print in this order
        reobsmaf.to_csv(site[k]+'_wtmp_ma_obs.csv',index=False,header=False,na_rep='NaN',float_format='%10.2f')
#########creat obs dc###############
        obsnewindex=[]
        for j in range(len(reobsda)):    
                obsnewindex.append(reobsda['mean'].index[j].replace(year=2000)) # puts all observations in the same year
        obsdc=pd.DataFrame(reobsda['mean'].values,index=obsnewindex)
#    tsdc=tsodnew.resample('D',how=['count','mean','median','min','max','std'],loffset=timedelta(hours=-12))
        reobsdc=obsdc[0].resample('D',how=['count','mean','median','min','max','std'])    #add columns for custom date format
        reobsdc['yy']=0
        reobsdc['mm']=reobsdc.index.month
        reobsdc['dd']=reobsdc.index.day
        output_fmt=['yy','mm','dd','count','mean','median','min','max','std']
        reobsdcf=reobsdc.reindex(columns=output_fmt)# found I needed to generate a new dataframe to print in this order
        reobsdcf.to_csv(site[k]+'_wtmp_dc_obs.csv',index=False,header=False,na_rep='NaN',float_format='%10.2f')
############creat obs mc###################
        obsmc=reobsdc['mean'].resample('m',how=['mean','median'],loffset=timedelta(days=-15))
        obsmc['count']=0
        obsmc['min']=0.
        obsmc['max']=0.
        obsmc['std']=0.
   
        ocount=reobsdc['count'].resample('m',how=['mean'],loffset=timedelta(days=-15)).values
        omi=reobsdc['min'].resample('m',how=['mean'],loffset=timedelta(days=-15)).values
        oma=reobsdc['max'].resample('m',how=['mean'],loffset=timedelta(days=-15)).values
        ostd=reobsdc['std'].resample('m',how=['mean'],loffset=timedelta(days=-15)).values
        for kk in range(len(obsmc)):
           obsmc['count'].values[kk]=ocount[kk]
           obsmc['min'].values[kk]=omi[kk]
           obsmc['max'].values[kk]=oma[kk]
           obsmc['std'].values[kk]=ostd[kk]
        obsmc['yy']=0
        obsmc['mm']=obsmc.index.month
        obsmc['dd']=0
        output_fmt=['yy','mm','dd','count','mean','median','min','max','std']
        obsmcf=obsmc.reindex(columns=output_fmt)# found I needed to generate a new dataframe to print in this order
        obsmcf.to_csv(site[k]+'_wtmp_mc_obs.csv',index=False,header=False,na_rep='NaN',float_format='%10.2f')

###################################################################

        urlfvcom = 'http://www.smast.umassd.edu:8080/thredds/dodsC/fvcom/hindcasts/30yr_gom3'
        nc = netCDF4.Dataset(urlfvcom)
        nc.variables
        lat = nc.variables['lat'][:]
        lon = nc.variables['lon'][:]
        times = nc.variables['time']
        jd = netCDF4.num2date(times[:],times.units)
        vname = 'temp'
        var = nc.variables[vname]

        inode = nearxy(lon,lat,loni,lati)
        modindex=netCDF4.date2index([obs_dt[0],obs_dt[-1]],times,select='nearest')
        modtso=pd.DataFrame(var[modindex[0]:modindex[1],44,inode],index=jd[modindex[0]:modindex[1]])
##############creat mod da#############
        remodda=float64(modtso[0]).resample('D',how=['count','mean','median','min','max','std'])
        remodda.ix[remodda['count']<minnumperday,['mean','median','min','max','std']] = 'NaN'
        remodda['yy']=remodda.index.year
        remodda['mm']=remodda.index.month
        remodda['dd']=remodda.index.day
        output_fmt=['yy','mm','dd','count','mean','median','min','max','std']
        remoddaf=remodda.reindex(columns=output_fmt)  
        remoddaf.to_csv(site[k]+'_wtmp_da_mod.csv',index=False,header=False,na_rep='NaN',float_format='%10.2f')
############creat mod ma#############
        remodma=float64(modtso[0]).resample('m',how=['count','mean','median','min','max','std'],kind='period')
        remodma.ix[remodma['count']<25*numperday,['mean','median','min','max','std']] = 'NaN'
        remodma['yy']=remodma.index.year
        remodma['mm']=remodma.index.month
        remodma['dd']=15
        output_fmt=['yy','mm','dd','count','mean','median','min','max','std']
        remodmaf=remodma.reindex(columns=output_fmt)# found I needed to generate a new dataframe to print in this order
        remodmaf.to_csv(site[k]+'_wtmp_ma_mod.csv',index=False,header=False,na_rep='NaN',float_format='%10.2f')
############creat mod dc##############
        newindex=[]
        for j in range(len(remodda)):    
                newindex.append(remodda['mean'].index[j].replace(year=2000)) # puts all observations in the same year
        moddc=pd.DataFrame(remodda['mean'].values,index=newindex)
#    tsdc=tsodnew.resample('D',how=['count','mean','median','min','max','std'],loffset=timedelta(hours=-12))
        remoddc=moddc[0].resample('D',how=['count','mean','median','min','max','std'])    #add columns for custom date format
        remoddc['yy']=0
        remoddc['mm']=remoddc.index.month
        remoddc['dd']=remoddc.index.day
        output_fmt=['yy','mm','dd','count','mean','median','min','max','std']
        remoddcf=remoddc.reindex(columns=output_fmt)# found I needed to generate a new dataframe to print in this order
        remoddcf.to_csv(site[k]+'_wtmp_dc_mod.csv',index=False,header=False,na_rep='NaN',float_format='%10.2f')
############creat mod mc##############
        modmc=remoddc['mean'].resample('m',how=['mean','median'],loffset=timedelta(days=-15))
        modmc['count']=0
        modmc['min']=0.
        modmc['max']=0.
        modmc['std']=0.
   
        tcount=remoddc['count'].resample('m',how=['mean'],loffset=timedelta(days=-15)).values
        tmi=remoddc['min'].resample('m',how=['mean'],loffset=timedelta(days=-15)).values
        tma=remoddc['max'].resample('m',how=['mean'],loffset=timedelta(days=-15)).values
        tstd=remoddc['std'].resample('m',how=['mean'],loffset=timedelta(days=-15)).values
        for kk in range(len(modmc)):
           modmc['count'].values[kk]=tcount[kk]
           modmc['min'].values[kk]=tmi[kk]
           modmc['max'].values[kk]=tma[kk]
           modmc['std'].values[kk]=tstd[kk]
        modmc['yy']=0
        modmc['mm']=modmc.index.month
        modmc['dd']=0
        output_fmt=['yy','mm','dd','count','mean','median','min','max','std']
        modmcf=modmc.reindex(columns=output_fmt)# found I needed to generate a new dataframe to print in this order
        modmcf.to_csv(site[k]+'_wtmp_mc_mod.csv',index=False,header=False,na_rep='NaN',float_format='%10.2f')
##############plot da compare figure##################
        fig=plt.figure(figsize=(16,10))
        ax=fig.add_subplot(211)
        ax.plot_date(obs_dt,obs_temp,fmt='-')
        plt.grid()
        ax.plot_date(jd[modindex[0]:modindex[1]],var[modindex[0]:modindex[1],44,inode],fmt='-',color='red')#bottom most value equals 44
        plt.ylabel(var.units)
        plt.title('eMOLT site '+site[k]+' temp vs FVCOM '+'%s at node=%d (Lon=%.4f, Lat=%.4f)' % (vname, inode+1, lon[inode], lat[inode]))
        plt.legend(['observed','modeled'],loc='best')
        plt.show()
###############plot mc compare figure##################
        TimeDelta=obsmcf.index[-1]-obsmcf.index[0]          
        ax1 = fig.add_subplot(212)
        my_x_axis_format(ax1, TimeDelta)
        ax1.plot_date(obsmcf.index,obsmcf['mean'],fmt='-')
        plt.grid()
        ax1.plot_date(modmcf.index,modmcf['mean'],fmt='-',color='red')#bottom most value equals 44
        plt.ylabel(var.units)
        plt.title('eMOLT site '+site[k]+' temp vs FVCOM '+'%s at node=%d (Lon=%.4f, Lat=%.4f)' % (vname, inode+1, lon[inode], lat[inode]))
        plt.legend(['observed','modeled'],loc='best')
        plt.show()
        plt.savefig(site[k]+'da_mc_mod_obs.png')
############calculate the different#######################
        output_fmt=[0,'mean','median','min','max','std']
        diffmc=obsmcf-modmcf
        month=pd.DataFrame(range(1,13),index=diffmc.index)
        diffmc=diffmc.join(month)
        diffmcf=diffmc.reindex(columns=output_fmt)
        diffmcf.columns=['month','mean','median','min','max','std']
        diffmcf.to_csv(site[k]+'_wtmp_mc_mod_mc_obs.csv',index=False,header=True,na_rep='NaN',float_format='%10.2f')

        diffma=reobsmaf-remodmaf
        date=[]
        for i in range(len(diffma)):
            date.append(str(diffma.index[i].year)+'-'+str(diffma.index[i].month))
        datetimepd=pd.DataFrame(date,index=diffma.index)
        diffma=diffma.join(datetimepd)
        diffmaf=diffma.reindex(columns=output_fmt)
        diffmaf.columns=['Year-Month','mean','median','min','max','std']
        diffmaf.to_csv(site[k]+'_wtmp_ma_mod_ma_obs.csv',index=False,header=True,na_rep='NaN',float_format='%10.2f')

        diffda=reobsdaf-remoddaf
        day=[]
        for i in range(len(diffda)):
            day.append(str(diffda.index[i].year)+'-'+str(diffda.index[i].month)+'-'+str(diffda.index[i].day))
        daytime=pd.DataFrame(day,index=diffda.index)
        diffda=diffda.join(daytime)
        diffdaf=diffda.reindex(columns=output_fmt)
        diffdaf.columns=['date','mean','median','min','max','std']
        diffdaf.to_csv(site[k]+'_wtmp_da_mod_da_obs.csv',index=False,header=True,na_rep='NaN',float_format='%10.2f')

        diffdc=reobsdcf-remoddcf
        daydc=[]
        for i in range(len(diffdc)):
            daydc.append(str(diffdc.index[i].year)+'-'+str(diffdc.index[i].month)+'-'+str(diffdc.index[i].day))
        daydctime=pd.DataFrame(daydc,index=diffdc.index)
        diffdc=diffdc.join(daydctime)
        diffdcf=diffdc.reindex(columns=output_fmt)
        diffdcf.columns=['date','mean','median','min','max','std']
        diffdcf.to_csv(site[k]+'_wtmp_dc_mod_dc_obs.csv',index=False,header=True,na_rep='NaN',float_format='%10.2f')
        
#####################
        del lati,loni,on,obs_dt,obs_temp,obstso,reobsda,reobsdaf,reobsma,reobsmaf,obsdc,reobsdc,reobsdcf,obsmc,obsmcf,nc,lat,lon,times,jd,var,inode,
        modindex,modtso,remodda,remoddaf,moddc,remoddc,remoddcf,modmc,modmcf,fig,ax,diffmc,diffmcf,diffma,diffmaf,diffda,diffdaf,diffdc,diffdcf,datetime
        