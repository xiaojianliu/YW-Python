# -*- coding: utf-8 -*-
"""
Created on Thu Aug  8 13:52:06 2013

@author: jmanning
"""
'''
import numpy as np
import matplotlib.pyplot as plt
import itertools
r = 3
c = 4
x = np.linspace(0,c,c+1)
y = np.linspace(0,r,r+1)
pts = itertools.product(x,y)
plt.scatter(*zip(*pts), marker = 'o', s = 30, color = 'red')
X, Y = np.meshgrid(x, y)
deg = np.arctan(Y**3-3*Y-X)
QP = plt.quiver(X, Y, np.cos(deg), np.sin(deg))
plt.grid()
plt.show()
###################################
from itertools import product
import matplotlib.pyplot as plt
points = np.array(list(product(range(3),range(4))))
plt.plot(points[:,0],points[:,1],'ro')
plt.show()
'''
import pylab
import pandas as pd
import numpy as np
import datetime
import matplotlib.pyplot as plt
from mpl_toolkits.basemap import Basemap
import os
import math
def dist(lat1,lon1,lat2,lon2):
    radius = 6371 # km
    

    dlat = math.radians(lat2-lat1)
    dlon = math.radians(lon2-lon1)
    a = math.sin(dlat/2) * math.sin(dlat/2) + math.cos(math.radians(lat1)) \
        * math.cos(math.radians(lat2)) * math.sin(dlon/2) * math.sin(dlon/2)
    c = 2 * math.atan2(math.sqrt(a), math.sqrt(1-a))
    d = radius * c
    
    def calcBearing(lat1, lon1, lat2, lon2):
       dLon = lon2 - lon1
       y = math.sin(dLon) * math.cos(lat2)
       x = math.cos(lat1) * math.sin(lat2) \
           - math.sin(lat1) * math.cos(lat2) * math.cos(dLon)
       return math.atan2(y, x)
       
    bear= math.degrees(calcBearing(lat1, lon1, lat2, lon2))  
    return d,bear
'''
def dist(lat1, lon1, lat2, lon2):
    if 1000 > lon1 > 0:
        (lat1,lon1) = dd2dm(lat1,lon1)
    if 1000 > lon2 > 0: 
        (lat2,lon2) = dd2dm(lat2,lon2)
    pid180 = math.pi/180
    alat = (lat1)*pid180
    alon = (lon1)*pid180
    blat = (lat2)*pid180
    blon = (lon2)*pid180
  
    cix = math.cos(alat)*math.cos(alat)*math.cos(alon-blon)+math.sin(alat)*math.sin(alat)
    distkmx = 6371*math.tan(math.sqrt(abs(1-cix*cix)/cix))
    ciy = math.cos(alat)*math.cos(blat)*math.cos(alon-alon)+math.sin(alat)*math.sin(blat) 
    distkmy = 6371*math.tan(math.sqrt(abs(1-ciy*ciy)/ciy))

    if distkmx > 0.001:   
        bear = abs(math.atan(distkmy/distkmx))
        if blon > alon and blat > alat:
            bear = 90-bear*180/math.pi 
        elif blon > alon and blat <= alat:
            bear = 90+bear*180/math.pi
        elif blon < alon and blat <= alat:        
            bear = 270-bear*180/math.pi
        elif blon < alon and blat > alat:  
            bear = 270+bear*180/math.pi

    else:
        if blat >= alat:
            bear = 0
        else:
            bear = 180

    ci = math.cos(alat)*math.cos(blat)*math.cos(alon-blon)+math.sin(alat)*math.sin(blat)
    distkm = 6371*math.tan(math.sqrt(abs(1-ci*ci)/ci))
    return distkm, bear
'''
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
def closetime(dfreal,dfmod):
         for i in range(len(dfreal)):
             distance=[]
             dffinalmod=[]
             for k in range(len(dfmod)):
                 distance.append(dfreal['time'].values[i]-dfmod['time'].values[k])
                 if dfreal['time'].values[i]-dfmod['time'].values[k]==min(distance):
                       dffinalmod.append(dfmod.values[k])
         return dffinalmod
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


latstart=np.array([])
lonstart=np.array([])
mean_distance=np.array([])
mean_std=np.array([])

f = open('distance.csv', 'w')
FList = np.genfromtxt('drift_data/FList.csv',dtype=None,names=['FNs'],delimiter=',')
FNS=list(FList['FNs'])
#FNS=['ID_48202.csv']
for kk in range(len(FNS)):
     averagefinal=[]
     stdfinal=[]
     print FNS[kk]
     conclusionfinal=pd.DataFrame(columns=['time','distkm','modlat','modlon','reallat','reallon'])
################get the real drift data##################
     df=pd.read_csv('drift_data/'+FNS[kk],index_col=0,skiprows=9,sep=',',names=['ID','TimeRD','TIME_GMT','YRDAY0_GMT','LON_DD','LAT_DD','TEMP','DEPTH'])
#     print 'start time is '+df['TIME_GMT'].values[0],'end time is '+df['TIME_GMT'].values[-1]
     d={'lat':df['LAT_DD'],'lon':df['LON_DD'],'time':df['TimeRD']}
     dfreal=pd.DataFrame(d)
     dtstart = datetime.datetime.strptime(df['TIME_GMT'].values[0], "%Y-%m-%d")
     dtend=datetime.datetime.strptime(df['TIME_GMT'].values[-1], "%Y-%m-%d")
     try:     
       if dtstart<datetime.datetime(2010, 12, 31, 0, 0):
     ################generate a timelist###############
         delta = datetime.timedelta(days=1)
         timelist=[]
         while dtstart<=dtend:
#            print dtstart.strftime("%Y-%m-%d")
            timelist.append(dtstart.strftime("%Y-%m-%d"))
            dtstart+=delta
         timelist=[val for val in timelist if val in df['TIME_GMT'].values]
         if   len(timelist)%2<>0:
            timelistloop=len(timelist)/2
         else:
            timelistloop=(len(timelist)-1)/2
##################get the model data########################################################################################
         for i in range(timelistloop):
           dtstartsep=timelist[2*i]
           dtendsep=timelist[2*i+2]
           
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
           
           dfsep=df.ix[df['TIME_GMT']==dtstartsep]
           lond=np.array([dfsep['LON_DD'].values[0]])
           latd=np.array([dfsep['LAT_DD'].values[0]])
           #latstart=np.append(latstart,latd)
           #lonstart=np.append(lonstart,lond)
           dtstartsep = datetime.datetime.strptime(dtstartsep, "%Y-%m-%d")
           startyear=dtstartsep.year
           startmonth=dtstartsep.month
           startday=dtstartsep.day
           dtendsep=datetime.datetime.strptime(dtendsep, "%Y-%m-%d")
           endyear=dtendsep.year
           endmonth=dtendsep.month
           endday=dtendsep.day
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
                ti=datetime.datetime.fromordinal(int(tn))
                YEAR=str(ti.year)
                MO=str(ti.month).zfill(2)
                DA=str(ti.day).zfill(2)
                hr=(tn-int(tn))*24
                HR=str(int(np.round(hr))).zfill(2)            
                TS=YEAR+MO+DA+HR+'0000'
#                print TS    
                FNU=FN0+'GOM3_'+YEAR+'/'+'u0/'+TS+'_u0.npy'
                FNV=FN0+'GOM3_'+YEAR+'/'+'v0/'+TS+'_v0.npy'
                u=np.load(FNU).flatten()
                v=np.load(FNV).flatten() 
                lont[kt],latt[kt]=lond,latd
                lond,latd=RungeKutta4_lonlat(lond,latd,Grid,u,v,tau)
           toc=os.times()
#          
           d={'time':tt,'lat':latt,'lon':lont}
           dfmod=pd.DataFrame(d)
           dfrealfinal = dfreal[(dfreal.time < t1) & (dfreal.time >= t0)]
           dfmodfinal=[]
           dfmodfinal.append(dfmod.values[0])
           for i in range(1,len(dfrealfinal)):
               distance=[]
               for k in range(1,len(dfmod)):
                  distance.append(abs(dfrealfinal['time'].values[i]-dfmod['time'].values[k]))
               n=np.argmin(distance)
               dfmodfinal.append(dfmod.values[n])
           dfmodfinal=pd.DataFrame(dfmodfinal)
           distkmlist=[]
           for v in range(len(dfrealfinal['lat'].values)):
                  distkm, bear=dist(dfrealfinal['lat'].values[v], dfrealfinal['lon'].values[v], dfmodfinal[0].values[v], dfmodfinal[1].values[v])
                  #if math.isnan(distkm)==False:            
                  distkmlist.append(distkm)                  
                  print "the distance of two point is "+str(distkm)
           d={'reallat':dfrealfinal['lat'].values,'reallon':dfrealfinal['lon'].values,'modlat':dfmodfinal[0].values,'modlon':dfmodfinal[1].values,'distkm':distkmlist,'time':dfrealfinal['time'].values}
           conclusion=pd.DataFrame(d)
           conclusionfinal=conclusionfinal.append(conclusion)
           sums=conclusion.distkm.values.sum()
           stdivi=conclusion.distkm.std()
           average=sums/len(conclusion)
           if math.isnan(average)==False and math.isnan(stdivi)==False and average<=2*stdivi and average<=50: 
               mean_distance=np.append(mean_distance,average)
               mean_std=np.append(mean_std,stdivi)
               londstart=np.array([dfsep['LON_DD'].values[0]])
               latdstart=np.array([dfsep['LAT_DD'].values[0]])
               latstart=np.append(latstart,latdstart)
               lonstart=np.append(lonstart,londstart)
           stdivi=conclusion.distkm.std()
           averagefinal.append(average)
           stdfinal.append(stdivi)
     except:
         kk+=1
     
'''
    ix=np.digitize(x,xbins)
    iy=np.digitize(y,ybins)
    xb=0.5*(xbins[:-1]+xbins[1:]) # bin x centers
    yb=0.5*(ybins[:-1]+ybins[1:]) # bin y centers
    zb_mean=np.empty((len(xbins)-1,len(ybins)-1),dtype=z.dtype)
    zb_median=np.empty((len(xbins)-1,len(ybins)-1),dtype=z.dtype)
    zb_std=np.empty((len(xbins)-1,len(ybins)-1),dtype=z.dtype)
    zb_num=np.zeros((len(xbins)-1,len(ybins)-1),dtype=int)    
    for iix in range(1,len(xbins)):
        for iiy in range(1,len(ybins)):
            k,=np.where((ix==iix) & (iy==iiy))
            zb_mean[iix-1,iiy-1]=np.mean(z[k])
            zb_median[iix-1,iiy-1]=np.median(z[k])
            zb_std[iix-1,iiy-1]=np.std(z[k])
            zb_num[iix-1,iiy-1]=len(z[k])
'''
xi = np.arange(-72,-66,0.25)
yi = np.arange(40.,45.,0.25)
plt.figure(figsize=(7,6))
#########bathymetry###############
latmap=[40.0,45.0]
lonmap=[-71.,-65.0]
url='http://geoport.whoi.edu/thredds/dodsC/bathy/gom03_v1_0'
from pydap.client import open_url
dataset=open_url(url)
basemap_lat=dataset['lat']
basemap_lon=dataset['lon']
basemap_topo=dataset['topo']

minlat=min(latmap)-0.01
maxlat=max(latmap)+0.01
minlon=min(lonmap)+0.01
maxlon=max(lonmap)-0.01
index_minlat=int(round(np.interp(minlat,basemap_lat,range(0,basemap_lat.shape[0]))))
index_maxlat=int(round(np.interp(maxlat,basemap_lat,range(0,basemap_lat.shape[0]))))
index_minlon=int(round(np.interp(minlon,basemap_lon,range(0,basemap_lon.shape[0]))))
index_maxlon=int(round(np.interp(maxlon,basemap_lon,range(0,basemap_lon.shape[0]))))
min_index_lat=min(index_minlat,index_maxlat)
max_index_lat=max(index_minlat,index_maxlat)
min_index_lon=min(index_minlon,index_maxlon)
max_index_lon=max(index_minlon,index_maxlon)

if index_minlat==0 or index_maxlat==0 or index_minlon==0 or index_maxlon==0:
        url='http://geoport.whoi.edu/thredds/dodsC/bathy/crm_vol1.nc'
        dataset=open_url(url)
        basemap_lat=dataset['lat']
        basemap_lon=dataset['lon']
        basemap_topo=dataset['topo']
        minlat=min(latmap)-0.01
        maxlat=max(latmap)+0.01
        minlon=min(lonmap)+0.01
        maxlon=max(lonmap)-0.01
        basemap_lat=[float(i) for i in basemap_lat]
        basemap_lat.reverse()
        range_basemap_lat=range(len(basemap_lat))
        range_basemap_lat.reverse()
        index_minlat=int(round(np.interp(minlat,basemap_lat,range_basemap_lat)))
        index_maxlat=int(round(np.interp(maxlat,basemap_lat,range_basemap_lat)))
        index_minlon=int(round(np.interp(minlon,basemap_lon,range(0,basemap_lon.shape[0]))))
        index_maxlon=int(round(np.interp(maxlon,basemap_lon,range(0,basemap_lon.shape[0]))))
min_index_lat=min(index_minlat,index_maxlat)
max_index_lat=max(index_minlat,index_maxlat)
X,Y=np.meshgrid(basemap_lon[index_minlon:index_maxlon],basemap_lat[min_index_lat:max_index_lat])
CS=plt.contour(X,Y,basemap_topo.topo[min_index_lat:max_index_lat,index_minlon:index_maxlon],[-70,-150],colors='black',linewith=0.2)
plt.clabel(CS, fontsize=7,fmt='%5.0f', inline=1)
#################geographic map##################################
#plt.figure(figsize=(7,6))

m = Basemap(projection='cyl',llcrnrlat=40.125,urcrnrlat=44.625,\
            llcrnrlon=-71.875,urcrnrlon=-66.375,resolution='h')
m.drawparallels(np.arange(int(40.125),int(44.625),1),labels=[1,0,0,0],linewidth=0.0)
m.drawmeridians(np.arange(int(-71.875),int(-66.375),1),labels=[0,0,0,1],linewidth=0.0)
m.drawcoastlines()
m.fillcontinents(color='grey')
m.drawmapboundary()
ix=np.digitize(lonstart,xi)
iy=np.digitize(latstart,yi)
zb_mean=np.empty((len(xi),len(yi)),dtype=mean_distance.dtype)
zb_count=np.empty((len(xi),len(yi)),dtype=mean_distance.dtype)
for iix in range(1,len(xi)):
    for iiy in range(1,len(yi)):
        k,=np.where((ix==iix) & (iy==iiy))
        zb_mean[iix-1,iiy-1]=np.mean(mean_distance[k-1])#you cant use k,if you do will rise warning like "index 12206 out of bounds 0<=index<12206"                
        zb_count[iix-1,iiy-1]=len(k)
#zb_mean=np.ma.masked_invalid(zb_mean)
xb=0.5*(xi[:-1]+xi[1:]) # bin x centers
yb=0.5*(yi[:-1]+yi[1:]) # bin y centers
xxb,yyb = np.meshgrid(xb, yb)
#im = NonUniformImage(ax, interpolation='nearest', extent=(-71.875,-66.375,40.625,45.125),
 #                    cmap=cm.Purples)
#im.set_data(xb, yb,zb_mean.T[:-1,:-1])

#ax.images.append(im)
#c=plt.pcolor(xxb,yyb,zb_mean.T,vmin=0,vmax=30,interpolation='nearest')
#plt.title('separate_distance',fontsize=18)
#plt.clim(0,35)
#b = plt.colorbar(c, orientation='vertical')
#plt.tick_params(axis='both', which='major', labelsize=20)
#plt.show() 
''' 
from scipy import interpolate
from scipy.interpolate import Rbf
from matplotlib import cm

xnew,ynew = np.mgrid[-71:-65.0:70j,40.0:45.0:70j]
tck = interpolate.bisplrep(xxb,yyb,zb_mean.T[:-1,:-1],s=0)
zb_mean_new = interpolate.bisplev(xnew[:,0],ynew[0,:],tck)
plt.figure()
plt.pcolor(xxb,yyb,zb_mean_new,vmin=0,vmax=35)
plt.colorbar()



'''

levels = [0, 3, 6, 9, 12, 15, 18, 21,24,27,30,33]
Q=plt.contourf(xxb,yyb,zb_mean.T[:-1,:-1],levels=levels,interpolation='bilinear')
thismanager4 = pylab.get_current_fig_manager()
plt.title('separate_distance',fontsize=18)
b = plt.colorbar(Q, orientation='vertical')
#b=mpl.colorbar.ColorbarBase(Q,norm=mpl.colors.Normalize(vmin=0, vmax=50))
plt.show()    
'''
fig = plt.figure()
fig.suptitle('separate_distance')
ax = fig.add_subplot(111)
from matplotlib.image import NonUniformImage
from matplotlib import cm
im = NonUniformImage(ax, interpolation='bilinear', extent=(-71.,-65.0,40.0,45.0),
                     cmap=cm.gray,vmin=0, vmax=40)
im.set_data(xb, yb, zb_mean.T[:-1,:-1])
ax.images.append(im)
ax.set_xlim(-71.,-65.0)
ax.set_ylim(40.0,45.0)
plt.colorbar()
plt.show()
'''


########number grid#################33
plt.figure(figsize=(7,6))
m = Basemap(projection='cyl',llcrnrlat=40,urcrnrlat=45,\
            llcrnrlon=-72,urcrnrlon=-66,resolution='h')
m.drawparallels(np.arange(int(40),int(45),1),labels=[1,0,0,0])
m.drawmeridians(np.arange(int(-72),int(-66),1),labels=[0,0,0,1])
m.drawparallels(np.arange(int(40),int(45),0.25))#,labels=[1,0,0,0])
m.drawmeridians(np.arange(int(-72),int(-66),0.25))

'''
m = Basemap(projection='cyl',llcrnrlat=40.5,urcrnrlat=45.5,\
            llcrnrlon=-72,urcrnrlon=-66,resolution='h')
m.drawparallels(np.arange(int(40.5),int(46.5),0.25))#,labels=[1,0,0,0])
m.drawmeridians(np.arange(int(-72),int(-66),0.25))#,labels=[0,0,0,1])
m.drawparallels(np.arange(int(40.5),int(45.5),1),labels=[1,0,0,0])
m.drawmeridians(np.arange(int(-72),int(-66),1),labels=[0,0,0,1])
'''
m.drawcoastlines()
m.fillcontinents(color='grey')
m.drawmapboundary()
for x in range(0,23):
    for y in range(0,19):
        if math.isnan(zb_mean[x,y])==False: 
             plt.text(xxb[y,x],yyb[y,x],round(zb_mean[x,y],1),fontsize=10,fontweight='bold',
                   ha='center',va='center',color='b')
             #plt.text(xxb[y,x],yyb[y,x],zb_count[x,y],fontsize=10,fontweight='bold',
                   #ha='center',va='center',color='b')
plt.title('separate_distance',fontsize=18)
plt.tick_params(axis='both', which='major', labelsize=20)
plt.show()