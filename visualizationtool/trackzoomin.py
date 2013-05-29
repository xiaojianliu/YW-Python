import matplotlib.pyplot as plt
from mpl_toolkits.basemap import Basemap
import numpy as np
import pylab
from datetime import datetime
from pydap.client import open_url
from datetime import timedelta
import sys

#############get the index of lat and lon???
def nearlonlat(lon,lat,lonp,latp):
    cp=np.cos(latp*np.pi/180.)
    # approximation for small distance
    dx=(lon-lonp)*cp
    dy=lat-latp
    dist2=dx*dx+dy*dy
    #dist1=np.abs(dx)+np.abs(dy)
    i=np.argmin(dist2)
    #min_dist=np.sqrt(dist2[i])
    return i#,min_dist 
    

urlname=open("ctrl_trackzoomin.csv", "r").readlines()[0][27:-1]
depth=int(open("ctrl_trackzoomin.csv", "r").readlines()[1][22:-1])
TIME=open("ctrl_trackzoomin.csv", "r").readlines()[2][31:-1]
numdays=int(open("ctrl_trackzoomin.csv", "r").readlines()[3][24:-1])


if urlname=="massbay":
    TIME=datetime.strptime(TIME, "%Y-%m-%d %H:%M:%S") 
    now=datetime.now()
    if TIME>now:
         diff=(TIME-now).days
    else:
         diff=(now-TIME).days
    if diff>3:
        print "please check your input start time,within 3 days both side form now on"
        sys.exit(0)
    numdays=timedelta(days=numdays)
    if TIME+numdays>now+timedelta(days=3):
        print "please check your numday.access period is between [now-3days,now+3days]"
        sys.exit(0)

latsize=[39,45]
lonsize=[-72.,-66]
fig=plt.figure(figsize=(7,6))
m = Basemap(projection='cyl',llcrnrlat=min(latsize)-0.01,urcrnrlat=max(latsize)+0.01,\
            llcrnrlon=min(lonsize)-0.01,urcrnrlon=max(lonsize)+0.01,resolution='h')#,fix_aspect=False)
m.drawparallels(np.arange(int(min(latsize)),int(max(latsize))+1,1),labels=[1,0,0,0])
m.drawmeridians(np.arange(int(min(lonsize)),int(max(lonsize))+1,1),labels=[0,0,0,1])
m.drawcoastlines()
m.fillcontinents(color='grey')
m.drawmapboundary()
if urlname=='30yr':
     stime=datetime.strptime(TIME, "%Y-%m-%d %H:%M:%S")
     timesnum=stime.year-1981
     standardtime=datetime.strptime(str(stime.year)+'-01-01 00:00:00', "%Y-%m-%d %H:%M:%S")
     timedeltaprocess=(stime-standardtime).days
     startrecord=26340+35112*(timesnum/4)+8772*(timesnum%4)+1+timedeltaprocess*24
     endrecord=startrecord+24*numdays
     url='http://www.smast.umassd.edu:8080/thredds/dodsC/fvcom/hindcasts/30yr_gom3?'+'lon,lat,latc,lonc,siglay,h,Times['+str(startrecord)+':1:'+str(startrecord)+']'
else:
     timeperiod=(TIME+numdays)-(now-timedelta(days=3))
     startrecord=(timeperiod.seconds)/60/60
     endrecord=startrecord+24*(numdays.days)
     url="http://www.smast.umassd.edu:8080/thredds/dodsC/FVCOM/NECOFS/Forecasts/NECOFS_FVCOM_OCEAN_MASSBAY_FORECAST.nc?"+'lon,lat,latc,lonc,siglay,h,Times['+str(startrecord)+':1:'+str(startrecord)+']'
dataset = open_url(url)
latc = np.array(dataset['latc'])
lonc = np.array(dataset['lonc'])
lat = np.array(dataset['lat'])
lon = np.array(dataset['lon'])
siglay=np.array(dataset['siglay'])
h=np.array(dataset['h'])
m.plot(lon,lat,'r.',lonc,latc,'b+')

###############################################################################    
def onclick(event):
    print 'button=%d, x=%d, y=%d, xdata=%f, ydata=%f'%(
        event.button, event.x, event.y, event.xdata, event.ydata)
    if event.button==3:
        latsize=[event.ydata-0.6,event.ydata+0.6]
        lonsize=[event.xdata-0.6,event.xdata+0.6]
        plt.figure(figsize=(7,6))
        m = Basemap(projection='cyl',llcrnrlat=min(latsize)-0.01,urcrnrlat=max(latsize)+0.01,\
            llcrnrlon=min(lonsize)-0.01,urcrnrlon=max(lonsize)+0.01,resolution='h')#,fix_aspect=False)
        m.drawparallels(np.arange(int(min(latsize)),int(max(latsize))+1,1),labels=[1,0,0,0])
        m.drawmeridians(np.arange(int(min(lonsize)),int(max(lonsize))+1,1),labels=[0,0,0,1])
        m.drawcoastlines()
        m.fillcontinents(color='grey')
        m.drawmapboundary()
        m.plot(lon,lat,'r.',lonc,latc,'b+')
        plt.show()
        spoint = pylab.ginput(1)
        la=spoint[0][1]
        lo=spoint[0][0]
        latd,lond=[],[]

        kf=nearlonlat(lonc,latc,lo,la) # nearest triangle center F - face
        kv=nearlonlat(lon,lat,lo,la)
        depthtotal=siglay[:,kv]*h[kv]
        layer=np.argmin(abs(depthtotal-depth))
        for i in range(startrecord,endrecord):
############read the particular time model from website#########    
               timeurl='['+str(i)+':1:'+str(i)+']'
               uvposition=str([layer])+str([kf])
               if urlname=="30yr":
                       url='http://www.smast.umassd.edu:8080/thredds/dodsC/fvcom/hindcasts/30yr_gom3?'+'Times'+timeurl+',u'+timeurl+uvposition+','+'v'+timeurl+uvposition 
               else:
                       url="http://www.smast.umassd.edu:8080/thredds/dodsC/FVCOM/NECOFS/Forecasts/NECOFS_FVCOM_OCEAN_MASSBAY_FORECAST.nc?"+'Times'+timeurl+',u'+timeurl+uvposition+','+'v'+timeurl+uvposition
                
               dataset = open_url(url)
               u=np.array(dataset['u'])
               v=np.array(dataset['v'])  
################get the point according the position###################
               print kf,u[0,0,0],v[0,0,0],layer
               par_u=u[0,0,0]
               par_v=v[0,0,0]   
               xdelta=par_u*60*60
               ydelta=par_v*60*60
               latdelta=ydelta/111111
               londelta=(xdelta/(111111*np.cos(la*np.pi/180)))
               la=la+latdelta
               lo=lo+londelta
               latd.append(la)
               lond.append(lo)
               kf=nearlonlat(lonc,latc,lo,la) # nearest triangle center F - face
               kv=nearlonlat(lon,lat,lo,la)# nearest triangle vertex V - vertex 
               depthtotal=siglay[:,kv]*h[kv]
               layer=np.argmin(abs(depthtotal-depth))
        plt.figure(figsize=(7,6))
        latsize=[min(latd)-0.6,max(latd)+0.6]
        lonsize=[min(lond)-0.6,max(lond)+0.6]
        m = Basemap(projection='cyl',llcrnrlat=min(latsize)-0.01,urcrnrlat=max(latsize)+0.01,\
            llcrnrlon=min(lonsize)-0.01,urcrnrlon=max(lonsize)+0.01,resolution='h')#,fix_aspect=False)
        m.drawparallels(np.arange(int(min(latsize)),int(max(latsize))+1,1),labels=[1,0,0,0])
        m.drawmeridians(np.arange(int(min(lonsize)),int(max(lonsize))+1,1),labels=[0,0,0,1])
        m.drawcoastlines()
        m.fillcontinents(color='grey')
        m.drawmapboundary()
        m.plot(lon,lat,'r.',lonc,latc,'b+')
        m.plot(lond,latd,linewidth=1.5,color='r')
        plt.title(urlname+' model map track Depth:'+str(depth)+' Time:'+TIME) 
        plt.savefig(urlname+'driftrack.png')
    return True
cid= fig.canvas.mpl_connect('button_press_event', onclick)
plt.show()
