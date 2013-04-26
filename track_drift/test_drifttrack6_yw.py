# -*- coding: utf-8 -*-
"""
Drifter Tracking using velocity field from FVCOM GOM3 model
2013-04-11 ver1 Runge-Kutta scheme for 2D field
2013-04-12 ver2 xy coordinates

@author: Vitalii Sheremet, FATE Project, 2012-2013
"""

import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.basemap import Basemap
#from fate_triangle import *
#import os
#from SeaHorseLib import *
from datetime import *

def RataDie(Y,M,D):
    RD=367*Y-(7*(Y+((M+9)/12))/4)-(3*(((Y+(M-9)/7)/100)+1)/4)+(275*M/9)+D-396;
    return RD

def RungeKutta4_2D(x,y,Grid,u,v,dt):
    """
Use classical 4th order 4-stage Runge-Kutta algorithm 
to track particles one time step

    x,y=RungeKutta4_2D(x,y,Grid,u,v,dt)

    x,y coordinates of an array of particles
    Grid - triangular grid info
    u,v  - velocity field defined on the grid
    Veli - velocity field interpolating function
           u,v=Veli(x,y,Grid,u,v)

Vitalii Sheremet, FATE Project, 2012-2013
    """
        
    x1=x*1.;       y1=y*1.;       u1,v1=Veli(x1,y1,Grid,u,v);  
    x2=x+0.5*dt*u1;y2=y+0.5*dt*v1;u2,v2=Veli(x2,y2,Grid,u,v);
    x3=x+0.5*dt*u2;y3=y+0.5*dt*v2;u3,v3=Veli(x3,y3,Grid,u,v);
    x4=x+0.5*dt*u3;y4=y+0.5*dt*v3;u4,v4=Veli(x4,y4,Grid,u,v);
    x=x+dt/6.*(u1+2.*u2+2.*u3+u4);
    y=y+dt/6.*(v1+2.*v2+2.*v3+v4);
    return x,y

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
    #dist1=np.abs(dx)+np.abs(dy)
    i=np.argmin(dist2)
    #min_dist=np.sqrt(dist2[i])
    return i#,min_dist 
    
def nearxy(x,y,xp,yp):
    """
i=nearxy(x,y,xp,yp)
find the closest node in the array (x,y) to a point (xp,yp)
input:
x,y - np.arrays of the grid nodes, cartesian coordinates
xp,yp - point on a plane
output:
i - index of the closest node
min_dist - the distance to the closest node
For coordinates on a sphere use function nearlonlat

Vitalii Sheremet, FATE Project
    """
    dx=x-xp
    dy=y-yp
    dist2=dx*dx+dy*dy
# dist1=np.abs(dx)+np.abs(dy)
    i=np.argmin(dist2)
    return i


def find_kf(Grid,xp,yp):
    """
kf,lamb0,lamb1,lamb2=find_kf(Grid,xp,yp)

find to which triangle a point (xp,yp) belongs
input:
Grid - triangular grid info
xp,yp - point on a plane
output:
kf - index of the the triangle
lamb0,lamb1,lamb2 - barycentric coordinates of P in the triangle

Vitalii Sheremet, FATE Project
    """

# coordinates of the vertices
    kvf=Grid['kvf']
    x=Grid['x'][kvf];y=Grid['y'][kvf]  
# calculate baricentric trilinear coordinates
    A012=((x[1,:]-x[0,:])*(y[2,:]-y[0,:])-(x[2,:]-x[0,:])*(y[1,:]-y[0,:])) 
# A012 is twice the area of the whole triangle,
# or the determinant of the linear system above.
# When xc,yc is the baricenter, the three terms in the sum are equal.
# Note the cyclic permutation of the indices
    lamb0=((x[1,:]-xp)*(y[2,:]-yp)-(x[2,:]-xp)*(y[1,:]-yp))/A012
    lamb1=((x[2,:]-xp)*(y[0,:]-yp)-(x[0,:]-xp)*(y[2,:]-yp))/A012
    lamb2=((x[0,:]-xp)*(y[1,:]-yp)-(x[1,:]-xp)*(y[0,:]-yp))/A012
    kf,=np.argwhere((lamb0>=0.)*(lamb1>=0.)*(lamb2>=0.))
#    kf=np.argwhere((lamb0>=0.)*(lamb1>=0.)*(lamb2>=0.)).flatten()
#    kf,=np.where((lamb0>=0.)*(lamb1>=0.)*(lamb2>=0.))
    return kf,lamb0[kf],lamb1[kf],lamb2[kf]

def polygonal_barycentric_coordinates(xp,yp,xv,yv):
    """
Calculate generalized barycentric coordinates within an N-sided polygon.

    w=polygonal_barycentric_coordinates(xp,yp,xv,yv)
    
    xp,yp - a point within an N-sided polygon
    xv,yv - vertices of the N-sided polygon, length N
    w     - polygonal baricentric coordinates, length N,
            normalized w.sum()=1
   
Used for function interpolation:
    fp=(fv*w).sum()
    where fv - function values at vertices,
    fp the interpolated function at the point (xp,yp)
    
Vitalii Sheremet, FATE Project    
    """
    N=len(xv)   
    j=np.arange(N)
    ja=(j+1)%N # next vertex in the sequence 
    jb=(j-1)%N # previous vertex in the sequence
# area of the chord triangle j-1,j,j+1
    Ajab=np.cross(np.array([xv[ja]-xv[j],yv[ja]-yv[j]]).T,np.array([xv[jb]-xv[j],yv[jb]-yv[j]]).T) 
# area of triangle p,j,j+1
    Aj=np.cross(np.array([xv[j]-xp,yv[j]-yp]).T,np.array([xv[ja]-xp,yv[ja]-yp]).T)  

# In FVCOM A is O(1.e7 m2) .prod() may result in inf
# to avoid this scale A
    Aj=Aj/max(abs(Aj))
    Ajab=Ajab/max(abs(Ajab))
    
    w=xv*0.
    j2=np.arange(N-2)
    for j in range(N):
# (j2+j+1)%N - list of triangles except the two adjacent to the edge pj
# For hexagon N=6 j2=0,1,2,3; if j=3  (j2+j+1)%N=4,5,0,1
        w[j]=Ajab[j]*Aj[(j2+j+1)%N].prod()

# normalize w so that sum(w)=1       
    w=w/w.sum() 
       
    return w
    
def Veli(xp,yp,Grid,u,v):
    """
Velocity interpolatin function

    ui,vi=Veli(x,y,Grid,u,v)
    
    xp,yp - arrays of points where the interpolated velocity is desired
    Grid - parameters of the triangular grid
    u,v - velocity field defined at the triangle baricenters
    
    """
    
# find the nearest vertex    
    kv=nearxy(Grid['x'],Grid['y'],xp,yp)
#    print kv
# list of triangles surrounding the vertex kv    
    kfv=Grid['kfv'][0:Grid['nfv'][kv],kv]
#    print kfv
    xv=Grid['xc'][kfv];yv=Grid['yc'][kfv]
    w=polygonal_barycentric_coordinates(xp,yp,xv,yv)
#    print w

# interpolation within polygon, w - normalized weights: w.sum()=1.    
    ui=(u[kfv]*w).sum()
    vi=(v[kfv]*w).sum()
        
    return ui,vi
    
def InterpF2V(uf,vf,Grid):
    """
Velocity interpolatin function

    uv,vv=InterpF2V(uf,vf,Grid)
    
    uf,vf - arrays of velocity defined at triangle baricenter
    Grid - parameters of the triangular grid
    uv,vv - arrays of velocity defined at triangle vertices

    kfv - array(NFVMAX,NV) indices of Faces surrounding a Vertex
    nfv - array(NV,) number of Faces surrounding a Vertex
    
    """
#GOM3Grid['kfv'][:,2100]
#Out[40]: array([3824, 3825, 3826, 4310, 4309, 4308, 3824,   -1,   -1], dtype=int32)   
#GOM3Grid['nfv'][2100]
#Out[41]: 6

#kfv[0:nfv[2100],2100]
#Out[74]: array([3824, 3825, 3826, 4310, 4309, 4308], dtype=int32)
    
    kfv=Grid['kfv']
    nfv=Grid['nfv']

    uv=x*0.;vv=y*0.
    for kv in range(len(x)):
#        kfv[0:nfv[kv],kv] six triangles surrounding the node kv
        uv[kv]=uf[kfv[0:nfv[kv],kv]].mean()
        vv[kv]=vf[kfv[0:nfv[kv],kv]].mean()
    
    
    return uv,vv    

############################################################################
# Drifter Tracking Script
############################################################################
    
# FVCOM OM3 triangular grid
"""
from pydap.client import open_url
from netCDF4 import Dataset
URL='http://www.smast.umassd.edu:8080/thredds/dodsC/fvcom/hindcasts/30yr_gom3'
#http://www.smast.umassd.edu:8080/thredds/dodsC/fvcom/hindcasts/30yr_gom3?
#a1u[0:1:3][0:1:90414],a2u[0:1:3][0:1:90414],art1[0:1:48450],art2[0:1:48450],
#aw0[0:1:2][0:1:90414],awx[0:1:2][0:1:90414],awy[0:1:2][0:1:90414],cc_hvc[0:1:90414],
#h[0:1:48450],lat[0:1:48450],latc[0:1:90414],lon[0:1:48450],lonc[0:1:90414],
#nbe[0:1:2][0:1:90414],nbsn[0:1:10][0:1:48450],nbve[0:1:8][0:1:48450],
#nn_hvc[0:1:48450],nprocs,ntsn[0:1:48450],ntve[0:1:48450],nv[0:1:2][0:1:90414],
#partition[0:1:90414],siglay[0:1:44][0:1:48450],siglev[0:1:45][0:1:48450],
#x[0:1:48450],xc[0:1:90414],y[0:1:48450],yc[0:1:90414],z0b[0:1:90414],
#Itime[0:1:171882],Itime2[0:1:171882],Times[0:1:171882],file_date[0:1:171882],
#iint[0:1:171882],kh[0:1:171882][0:1:45][0:1:48450],
#km[0:1:171882][0:1:45][0:1:48450],kq[0:1:171882][0:1:45][0:1:48450],
#l[0:1:171882][0:1:45][0:1:48450],net_heat_flux[0:1:171882][0:1:48450],
#omega[0:1:171882][0:1:45][0:1:48450],q2[0:1:171882][0:1:45][0:1:48450],
#q2l[0:1:171882][0:1:45][0:1:48450],salinity[0:1:171882][0:1:44][0:1:48450],
#short_wave[0:1:171882][0:1:48450],temp[0:1:171882][0:1:44][0:1:48450],
#time[0:1:171882],u[0:1:171882][0:1:44][0:1:90414],ua[0:1:171882][0:1:90414],
#uwind_stress[0:1:171882][0:1:90414],v[0:1:171882][0:1:44][0:1:90414],
#va[0:1:171882][0:1:90414],vwind_stress[0:1:171882][0:1:90414],
#ww[0:1:171882][0:1:44][0:1:90414],zeta[0:1:171882][0:1:48450]

#ds=open_url(URL)                 # pydap version 
ds = Dataset(URL,'r').variables   # netCDF4 version

#xxx=ds['xxx']; np.save('gom3.xxx.npy',np.array(xxx))
a1u=ds['a1u']; np.save('gom3.a1u.npy',np.array(a1u))
a2u=ds['a2u']; np.save('gom3.a2u.npy',np.array(a2u))
art1=ds['art1']; np.save('gom3.art1.npy',np.array(art1))
art2=ds['art2']; np.save('gom3.art2.npy',np.array(art2))
aw0=ds['aw0']; np.save('gom3.aw0.npy',np.array(aw0))
awx=ds['awx']; np.save('gom3.awx.npy',np.array(awx))
awy=ds['awy']; np.save('gom3.awy.npy',np.array(awy))
cc_hvc=ds['cc_hvc']; np.save('gom3.cc_hvc.npy',np.array(cc_hvc))
    
h=ds['h']; np.save('gom3.h.npy',np.array(h))

lat=ds['lat']; np.save('gom3.lat.npy',np.array(lat))
lon=ds['lon']; np.save('gom3.lon.npy',np.array(lon))
latc=ds['latc']; np.save('gom3.latc.npy',np.array(latc))
lonc=ds['lonc']; np.save('gom3.lonc.npy',np.array(lonc))

nbe=ds['nbe']; np.save('gom3.nbe.npy',np.array(nbe))
nbsn=ds['nbsn']; np.save('gom3.nbsn.npy',np.array(nbsn))
nbve=ds['nbve']; np.save('gom3.nbve.npy',np.array(nbve))
nn_hvc=ds['nn_hvc']; np.save('gom3.nn_hvc.npy',np.array(nn_hvc))
nprocs=ds['nprocs']; np.save('gom3.nprocs.npy',np.array(nprocs))
ntsn=ds['ntsn']; np.save('gom3.ntsn.npy',np.array(ntsn))
ntve=ds['ntve']; np.save('gom3.ntve.npy',np.array(ntve))
nv=ds['nv']; np.save('gom3.nv.npy',np.array(nv))
partition=ds['partition']; np.save('gom3.partition.npy',np.array(partition))
siglay=ds['siglay']; np.save('gom3.siglay.npy',np.array(siglay))
siglev=ds['siglev']; np.save('gom3.siglev.npy',np.array(siglev))

x=ds['x']; np.save('gom3.x.npy',np.array(x))
xc=ds['xc']; np.save('gom3.xc.npy',np.array(xc))
y=ds['y']; np.save('gom3.y.npy',np.array(y))
yc=ds['yc']; np.save('gom3.yc.npy',np.array(yc))
"""    
    
x=np.load('gom3.x.npy')
y=np.load('gom3.y.npy')
xc=np.load('gom3.xc.npy')
yc=np.load('gom3.yc.npy')

lon=np.load('gom3.lon.npy')
lat=np.load('gom3.lat.npy')
lonc=np.load('gom3.lonc.npy')
latc=np.load('gom3.latc.npy')



#nv: Array of 32 bit Integers [three = 0..2][nele = 0..90414] 
#long_name: nodes surrounding element
#standard_name: face_node_connectivity
#start_index: 1
nv=np.load('gom3.nv.npy')
nv-=1 # convert from FORTRAN to python 0-based indexing
#kvf=nv

#nbe: Array of 32 bit Integers [three = 0..2][nele = 0..90414] 
# long_name: elements surrounding each element
nbe=np.load('gom3.nbe.npy')
nbe-=1 # convert from FORTRAN to python 0-based indexing
#kff=nbe

#nbsn: Array of 32 bit Integers [maxnode = 0..10][node = 0..48450]
#long_name: nodes surrounding each node
 # list of nodes surrounding a given node, 1st and last entries identical to make a closed loop
nbsn=np.load('gom3.nbsn.npy')
nbsn-=1 # convert from FORTRAN to python 0-based indexing
#kvv=nbsn

#ntsn: Array of 32 bit Integers [node = 0..48450]
#long_name: #nodes surrounding each node
 # the number of nodes surrounding a given node + 1, because 1st and last entries identical to make a closed loop
ntsn=np.load('gom3.ntsn.npy')
#nvv=ntsn

#nbve: Array of 32 bit Integers [maxelem = 0..8][node = 0..48450] 
#long_name: elems surrounding each node
# list of elements surrounding a given node, 1st and last entries identical to make a closed loop
nbve=np.load('gom3.nbve.npy')
nbve-=1 # convert from FORTRAN to python 0-based indexing
#kfv=nbve

#ntve: Array of 32 bit Integers [node = 0..48450] 
#long_name: #elems surrounding each node
# the number of elements surrounding a given node + 1, because 1st and last entries identical to make a closed loop
ntve=np.load('gom3.ntve.npy')
#nfv=ntve

Grid={'x':x,'y':y,'xc':xc,'yc':yc,'lon':lon,'lat':lat,'lonc':lonc,'latc':latc,'kvf':nv,'kff':nbe,'kvv':nbsn,'nvv':ntsn,'kfv':nbve,'nfv':ntve}


#########################################################################################
latsize=[39.0,45.0]
lonsize=[-72.,-66.0]
plt.figure(figsize=(7,6))
m = Basemap(projection='cyl',llcrnrlat=min(latsize)-0.01,urcrnrlat=max(latsize)+0.01,\
            llcrnrlon=min(lonsize)-0.01,urcrnrlon=max(lonsize)+0.01,resolution='h')#,fix_aspect=False)
m.drawparallels(np.arange(int(min(latsize)),int(max(latsize))+1,1),labels=[1,0,0,0])
m.drawmeridians(np.arange(int(min(lonsize)),int(max(lonsize))+1,1),labels=[0,0,0,1])
m.drawcoastlines()
m.fillcontinents(color='grey')
m.drawmapboundary()

#strat point
'''
spoint = pylab.ginput(1)
la=spoint[0][1]
lo=spoint[0][0]
inode=nearlonlat(lonc,latc,lo,la) 
xd=xc(inode)
yd=yc(inode)
'''
#########################################################################################

#u=np.load('20010101000000_ua.npy').flatten()
#v=np.load('20010101000000_va.npy').flatten()
FN0='/home/vsheremet/FATE/GOM3_DATA/'
FN0='/media/FreeAgent GoFlex Drive/GOM3_DATA/'


xd=np.array([1200000.])
yd=np.array([-200000.])

t0=RataDie(2005,4,1)
t1=RataDie(2005,4,31)
tt=np.arange(t0,t1,1./24.)
NT=len(tt)
xt=np.zeros(NT)
yt=np.zeros(NT)
dt=60*60.
latplot,lonplot=[],[] 
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
    
    xt[kt],yt[kt]=xd,yd
    xd,yd=RungeKutta4_2D(xd,yd,Grid,u,v,dt)
################################

# got track xt,yt
# convert to lont,latt
lont=xt*0.;latt=yt*0.;
for kt in range(NT):
    xd=xt[kt];yd=yt[kt];
# find the nearest vertex    
    kv=nearxy(Grid['x'],Grid['y'],xd,yd)
    xv=Grid['x'][kv];yv=Grid['y'][kv]
    lonv=Grid['lon'][kv];latv=Grid['lat'][kv]
#    lond=(xd-xv)/(111111.*np.cos(latv*np.pi/180.))+lonv
#    latd=(yd-yv)/111111.+latv
    latd=(yd-yv)/111111.+latv
    lond=(xd-xv)/(111111.*np.cos((yd-yv)/111111.*np.pi/180.))+lonv
    lont[kt]=lond
    latt[kt]=latd
    
    


 
"""
    ####################
    xdelta=dt/6.*(u1+2.*u2+2.*u3+u4)
    ydelta=dt/6.*(v1+2.*v2+2.*v3+v4)
    latdelta=ydelta/111111
    londelta=xdelta/(111111*cos(lat*np.pi/180))
    ####################
   
    la=la+latdelta
    lo=lo+londelta
    latplot.append(la)
    lonplot.append(lo)
"""    
###############################
m.plot(lont,latt,linewidth=1.5,color='r')
plt.show()
plt.figure(figsize=(7,6))
plt.plot(x,y,'r.',xc,yc,'b+');
plt.plot(xt,yt,'ko-')
'''

plt.figure()
plt.plot(lon,lat,'r.',lonc,latc,'b+');
plt.plot(lont,latt,'ko-')
'''
plt.show()    

