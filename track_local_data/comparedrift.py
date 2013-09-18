# -*- coding: utf-8 -*-
"""
Created on Thu May  2 14:36:18 2013

@author: jmanning
"""

import pandas as pd
import numpy as np
from datetime import datetime
import matplotlib.pyplot as plt
from mpl_toolkits.basemap import Basemap
import os
FList = np.genfromtxt('drift_data/Drift_list.csv',dtype=None,names=['FNs'],delimiter=',')
FNS=list(FList['FNs'])
for k in range(len(FNS)):
     print FNS[k]
     df=pd.read_csv('drift_data/'+FNS[k],index_col=0,skiprows=9,sep=',',names=['ID','TimeRD','TIME_GMT','YRDAY0_GMT','LON_DD','LAT_DD','TEMP','DEPTH'])
     print 'start time is '+df['TIME_GMT'].values[0],'end time is '+df['TIME_GMT'].values[-1]
     dtstart = datetime.strptime(df['TIME_GMT'].values[0], "%Y-%m-%d")
     dtend=datetime.strptime(df['TIME_GMT'].values[-1], "%Y-%m-%d")
     startlat=df['LAT_DD'].values[0]
     startlon=df['LON_DD'].values[0]
     print dtstart,dtend,startlat,startlon
     latmax=max(df['LAT_DD'].values)
     latmin=min(df['LAT_DD'].values)
     lonmax=max(df['LON_DD'].values)
     lonmin=min(df['LON_DD'].values)
     plt.figure(figsize=(7,6))
     m = Basemap(projection='cyl',llcrnrlat=latmin-2,urcrnrlat=latmax+2,\
           llcrnrlon=lonmin-2,urcrnrlon=lonmax+2,resolution='h')
     m.drawparallels(np.arange(int(latmin),int(latmax)+1,1),labels=[1,0,0,0])
     m.drawmeridians(np.arange(int(lonmin),int(lonmax)+1,2),labels=[0,0,0,2])
     m.drawcoastlines()
     m.fillcontinents(color='grey')
     m.drawmapboundary()
     x, y = m(df['LON_DD'].values,df['LAT_DD'].values)
     m.plot(x,y,linewidth=1.5,color='r')
     plt.title(FNS[k]+' mod vs obs')

#######################################################################################################################################################
#######################################################################################################################################################

     def RungeKutta4_lonlat(lon,lat,Grid,u,v,tau): 
         lon1=lon*1.;          lat1=lat*1.;        urc1,v1=VelInterp_lonlat(lon1,lat1,Grid,u,v);  
         lon2=lon+0.5*tau*urc1;lat2=lat+0.5*tau*v1;urc2,v2=VelInterp_lonlat(lon2,lat2,Grid,u,v);
         lon3=lon+0.5*tau*urc2;lat3=lat+0.5*tau*v2;urc3,v3=VelInterp_lonlat(lon3,lat3,Grid,u,v);
         lon4=lon+0.5*tau*urc3;lat4=lat+0.5*tau*v3;urc4,v4=VelInterp_lonlat(lon4,lat4,Grid,u,v);
         lon=lon+tau/6.*(urc1+2.*urc2+2.*urc3+urc4);
         lat=lat+tau/6.*(v1+2.*v2+2.*v3+v4);    
         return lon,lat
     def nearxy(x,y,xp,yp):
         dx=x-xp
         dy=y-yp
         dist2=dx*dx+dy*dy
         i=np.argmin(dist2)
         return i
     def nearlonlat(lon,lat,lonp,latp):
         cp=np.cos(latp*np.pi/180.)
         dx=(lon-lonp)*cp
         dy=lat-latp
         dist2=dx*dx+dy*dy
         i=np.argmin(dist2)
         return i
     def find_kf(Grid,xp,yp):
         kvf=Grid['kvf']
         x=Grid['x'][kvf];y=Grid['y'][kvf]  
         A012=((x[1,:]-x[0,:])*(y[2,:]-y[0,:])-(x[2,:]-x[0,:])*(y[1,:]-y[0,:])) 
         lamb0=((x[1,:]-xp)*(y[2,:]-yp)-(x[2,:]-xp)*(y[1,:]-yp))/A012
         lamb1=((x[2,:]-xp)*(y[0,:]-yp)-(x[0,:]-xp)*(y[2,:]-yp))/A012
         lamb2=((x[0,:]-xp)*(y[1,:]-yp)-(x[1,:]-xp)*(y[0,:]-yp))/A012
         kf,=np.argwhere((lamb0>=0.)*(lamb1>=0.)*(lamb2>=0.))
         return kf,lamb0[kf],lamb1[kf],lamb2[kf]
     def find_kf2(Grid,xp,yp): 
         kv=nearxy(Grid['x'],Grid['y'],xp,yp)
         kfv=Grid['kfv'][0:Grid['nfv'][kv],kv]
         kf=nearxy(Grid['xc'],Grid['yc'],xp,yp)
         kff=Grid['kff'][:,kf]
         kkf=np.concatenate((kfv,np.array([kf]),kff))
         kvf=Grid['kvf'][:,kkf]
         x=Grid['x'][kvf];y=Grid['y'][kvf]
         A012=((x[1,:]-x[0,:])*(y[2,:]-y[0,:])-(x[2,:]-x[0,:])*(y[1,:]-y[0,:])) 
         lamb0=((x[1,:]-xp)*(y[2,:]-yp)-(x[2,:]-xp)*(y[1,:]-yp))/A012
         lamb1=((x[2,:]-xp)*(y[0,:]-yp)-(x[0,:]-xp)*(y[2,:]-yp))/A012
         lamb2=((x[0,:]-xp)*(y[1,:]-yp)-(x[1,:]-xp)*(y[0,:]-yp))/A012
         kf=np.argwhere((lamb0>=0.)*(lamb1>=0.)*(lamb2>=0.)).flatten()
         kf=kf[0]
         return kkf[kf],lamb0[kf],lamb1[kf],lamb2[kf]
     def polygonal_barycentric_coordinates(xp,yp,xv,yv):
         N=len(xv)   
         j=np.arange(N)
         ja=(j+1)%N # next vertex in the sequence 
         jb=(j-1)%N # previous vertex in the sequence
         Ajab=np.cross(np.array([xv[ja]-xv[j],yv[ja]-yv[j]]).T,np.array([xv[jb]-xv[j],yv[jb]-yv[j]]).T) 
         Aj=np.cross(np.array([xv[j]-xp,yv[j]-yp]).T,np.array([xv[ja]-xp,yv[ja]-yp]).T)  
         Aj=Aj/max(abs(Aj))
         Ajab=Ajab/max(abs(Ajab))    
         w=xv*0.
         j2=np.arange(N-2)
         for j in range(N):
             w[j]=Ajab[j]*Aj[(j2+j+1)%N].prod()   
         w=w/w.sum()        
         return w,Aj
     def Veli(x,y,Grid,u,v):
         kf=nearxy(Grid['xc'],Grid['yc'],x,y)
         ui=u[kf]
         vi=v[kf]
         return ui,vi
     def Veli2(xp,yp,Grid,u,v):
         kv=nearxy(Grid['x'],Grid['y'],xp,yp)   
         kfv=Grid['kfv'][0:Grid['nfv'][kv],kv]
         xv=Grid['xc'][kfv];yv=Grid['yc'][kfv]
         w=polygonal_barycentric_coordinates(xp,yp,xv,yv)
         ui=(u[kfv]*w).sum()
         vi=(v[kfv]*w).sum()
         return ui,vi
         '''
     def VelInterp_lonlat(lonp,latp,Grid,u,v):
         kv=nearlonlat(Grid['lon'],Grid['lat'],lonp,latp)
         kfv=Grid['kfv'][0:Grid['nfv'][kv],kv]
         lonv=Grid['lonc'][kfv];latv=Grid['latc'][kfv] 
         w=polygonal_barycentric_coordinates(lonp,latp,lonv,latv)
         cv=Grid['coslatc'][kfv]    
         urci=(u[kfv]/cv*w).sum()
         vi=(v[kfv]*w).sum()
         return urci,vi
         '''
     def VelInterp_lonlat(lonp,latp,Grid,u,v):
          kv=nearlonlat(Grid['lon'],Grid['lat'],lonp,latp)
          kfv=Grid['kfv'][0:Grid['nfv'][kv],kv]
          lonv=Grid['lonc'][kfv];latv=Grid['latc'][kfv]
          w,Aj=polygonal_barycentric_coordinates(lonp,latp,lonv,latv)
          Aj=Aj/Aj.sum()
          if np.argwhere(Aj<0).flatten().size>0:
             for kv1 in Grid['kvv'][0:Grid['nvv'][kv],kv]:
                 kfv1=Grid['kfv'][0:Grid['nfv'][kv1],kv1]
                 lonv1=Grid['lonc'][kfv1];latv1=Grid['latc'][kfv1]
                 w1,Aj1=polygonal_barycentric_coordinates(lonp,latp,lonv1,latv1)
                 Aj1=Aj1/Aj1.sum()
                 if np.argwhere(Aj1<0).flatten().size==0:
                     w=w1;kfv=kfv1;kv=kv1;Aj=Aj1
          if np.argwhere(w<0).flatten().size>0:
             print kv,kfv,w
          cv=Grid['coslatc'][kfv]
          urci=(u[kfv]/cv*w).sum()
          vi=(v[kfv]*w).sum()
          return urci,vi 
     def RataDie(yr,mo=1,da=1,hr=0,mi=0,se=0):
         RD=367*yr-(7*(yr+((mo+9)//12))//4)-(3*(((yr+(mo-9)//7)//100)+1)//4)+(275*mo//9)+da-396+(hr*3600+mi*60+se)/86400.;
         return RD    
     x=np.load('gom3.x.npy')
     y=np.load('gom3.y.npy')
     xc=np.load('gom3.xc.npy')
     yc=np.load('gom3.yc.npy')
     lon=np.load('gom3.lon.npy')
     lat=np.load('gom3.lat.npy')
     lonc=np.load('gom3.lonc.npy')
     latc=np.load('gom3.latc.npy')

     coslat=np.cos(lat*np.pi/180.)
     coslatc=np.cos(latc*np.pi/180.)

     nv=np.load('gom3.nv.npy')
     nv-=1 # convert from FORTRAN to python 0-based indexing
     nbe=np.load('gom3.nbe.npy')
     nbe-=1 # convert from FORTRAN to python 0-based indexing
     nbsn=np.load('gom3.nbsn.npy')
     nbsn-=1 # convert from FORTRAN to python 0-based indexing
     ntsn=np.load('gom3.ntsn.npy')
     nbve=np.load('gom3.nbve.npy')
     nbve-=1 # convert from FORTRAN to python 0-based indexing
     ntve=np.load('gom3.ntve.npy')
     Grid={'x':x,'y':y,'xc':xc,'yc':yc,'lon':lon,'lat':lat,'lonc':lonc,'latc':latc,'coslat':coslat,'coslatc':coslatc,'kvf':nv,'kff':nbe,'kvv':nbsn,'nvv':ntsn,'kfv':nbve,'nfv':ntve}

     FN0='/media/FreeAgent GoFlex Drive/GOM3_DATA/'

     lond=np.array([startlon])
     latd=np.array([startlat])
     startyear=dtstart.year
     startmonth=dtstart.month
     startday=dtstart.day
     endyear=dtend.year
     endmonth=dtend.month
     endday=dtend.day
     t0=RataDie(startyear,startmonth,startday)
     t1=RataDie(endyear,endmonth,endday)
     tt=np.arange(t0,t1,1./24.)
     NT=len(tt)
     lont=np.zeros(NT)
     latt=np.zeros(NT)
     dt=60*60.
     tau=dt/111111. # deg per (velocityunits*dt)
     tic=os.times()
     for kt in range(NT):
          tRD=tt[kt]
          tn=np.round(tRD*24.)/24.
          ti=datetime.fromordinal(int(tn))
          YEAR=str(ti.year)
          MO=str(ti.month).zfill(2)
          DA=str(ti.day).zfill(2)
          hr=(tn-int(tn))*24
          HR=str(int(np.round(hr))).zfill(2)            
          TS=YEAR+MO+DA+HR+'0000'
          print TS    
          FNU=FN0+'GOM3_'+YEAR+'/'+'u0/'+TS+'_u0.npy'
          FNV=FN0+'GOM3_'+YEAR+'/'+'v0/'+TS+'_v0.npy'
          u=np.load(FNU).flatten()
          v=np.load(FNV).flatten() 
          lont[kt],latt[kt]=lond,latd
          lond,latd=RungeKutta4_lonlat(lond,latd,Grid,u,v,tau)
     toc=os.times()
     print 'timing [s] per step: ', (toc[0]-tic[0])/NT,(toc[1]-tic[1])/NT
     m.plot(lont,latt,linewidth=1.5,color='blue')
     plt.show()
     plt.savefig(FNS[k]+'mod_vs_obs.png')