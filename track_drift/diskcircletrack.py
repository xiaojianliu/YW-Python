import numpy as np
import matplotlib.pyplot as plt
import os
from datetime import *
from mpl_toolkits.basemap import Basemap

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
       
    return w

def VelInterp_lonlat(lonp,latp,Grid,u,v):
    kv=nearlonlat(Grid['lon'],Grid['lat'],lonp,latp)
    kfv=Grid['kfv'][0:Grid['nfv'][kv],kv]
    lonv=Grid['lonc'][kfv];latv=Grid['latc'][kfv] 
    w=polygonal_barycentric_coordinates(lonp,latp,lonv,latv)
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

lond=np.array([-67.28])
latd=np.array([41.35])
t0=RataDie(1985,6,18)
t1=RataDie(1985,7,18)
tt=np.arange(t0,t1,1./24.)
NT=len(tt)
lont=np.zeros(NT)
latt=np.zeros(NT)
dt=60*60.
tau=dt/111111.
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

plt.figure()

plt.plot(lon,lat,'r.',lonc,latc,'b+');#triangle model show
plt.plot(lont,latt,'ko-',lont[-1],latt[-1],'mo')#show the drift track

#################################

FN='XY.npz'
Z=np.load(FN)
x=Z['X']
y=Z['Y']
for nn in range(1,5):
    for n in range(0,20):
        lond=np.array([x[n,nn]])
        latd=np.array([y[n,nn]])
        print n,nn,lond,latd
        t0=RataDie(1985,6,18)
        t1=RataDie(1985,7,18)
        tt=np.arange(t0,t1,1./24.)
        NT=len(tt)
        lont=np.zeros(NT)
        latt=np.zeros(NT)
        dt=60*60.
        tau=dt/111111.
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
              print 'timing [s] per step: ', (toc[0]-tic[0])/NT,(toc[1]-tic[1])/NT,n,nn
        plt.plot(lon,lat,'r.',lonc,latc,'b+');#triangle model show
        plt.plot(lont,latt,'ko-',lont[-1],latt[-1],'mo')#show the drift track
        plt.show()    
