"""
9/JAN/2013
Description:comparing the drift and model track between 1988-1990 years
    
Import problem:
    
Response:
"""
#comparing the drift and model track between 1988-1990 years

#give a list of drifters,
#startdate is the first value which we calculate from drift.

#every numdays, give a new beginning 
# xiuling   4 Jan

import sys
sys.path.append("/home3/ocn/jmanning/py/jmanning/")
#sys.path.append("D:/py/jmanning/")

from models import *
from drifter import *
from getdata import *
from conversions import *
from matplotlib.dates import num2date,date2num
import matplotlib.pyplot as plt
from matplotlib.dates import  DateFormatter

#### HARD CODES#################################
days=6 #total number of days wanted to track
numdays=2 #every two days, plot one track for drift and model
hours=1 # hours to step through
depth=[0]#, -25,-50,-75]  in meters
#############################################

track_num=int(np.ceil(float(days)/numdays))# number of comparison's per drifter
daystep=float(hours)/24  #fraction of a day to step

driftid=[#'19899850']
         '19886981','19886970','19886971','19886973','19886974','19886975','19886976','19886977','19886978','19886979','19886982','19886983',
         '19886984','19886985','19886986','19886987','19886988','19886989','19899840','19899841','19899842','19899843','19899844','19899845', '19899846',
         '19899847','19899848','19899850','19899851','19899852','19899854']
  
#get the model
url_name="30yr"
url='http://www.smast.umassd.edu:8080/thredds/dodsC/fvcom/hindcasts/30yr_gom3'
dataset=open_url(url)
times=dataset['Times']
list_year=[1978, 1979, 1980, 1981, 1982, 1983, 1984, 1985, 1986, 1987, 1988, 1989, 1990,1991]#year list for the 30yr url
list_year_index=[0,8772+1,8772*2+1,8772*3+24+1,8772*4+24+1,8772*5+24+1,8772*6+24+1,8772*7+24*2+1,8772*8+24*2+1,8772*9+24*2+1,8772*10+24*2+1,8772*11+24*3+1,
                 8772*12+24*3+1,8772*13+24*3+1,times.shape[0]]#each year start index
#open the html, and write the code to html
f_open=open("/net/nwebserver/epd/ocean/MainPage/fate/index.html","a")


'''path changed because in different computer'''
#f_open=open("D:/py/index.html","a")
f_open.write('<html><center><h1>Drifter VS model</h1><table border="1">')
f_open.write('<b><tr><b><td align="right">drift id</td></b>'
        '<b><td>SP<br>(mean)</td></b><b><td>std</td></b><b><td>case</td></b><b><td>depth</td></b>'
        '</tr></b>')

fig3=plt.figure(1)#plot the figure for total mean of distance
ax3= fig3.add_subplot(111)

all_drift_std,all_drift_mean,each_hour_days_dis=[],[],[]
for f in driftid:
  d_lat,d_lon,time_drift_date,dep=getdrift(f)#use function getdrift, get the drift data.
  #convert time_drift_date to number values
  time_drift=[date2num(i) for i in time_drift_date]
  '''  print d_yearday[0]
  # if the d_yearday is matlab time, convert to python time
  if d_yearday[0].month<>dt.datetime.strptime(d_time[0],'%Y-%m-%d').month or d_yearday[0].day<>dt.datetime.strptime(d_time[0],'%Y-%m-%d').day:
      d_yearday=[float(date2num(i))+1 for i in d_yearday]

  time_drift=[]#convert yearday to datetime
  for i in range(len(d_yearday)):
     time_drift.append(date2num(num2date(d_yearday[i]).replace(year=int(d_time[i][0:4]))))'''
  
  #if list of d_yearday is not sorted, we need sort it.
  index_d_yearday=range(len(time_drift))
  index_d_yearday.sort(lambda x,y:cmp(time_drift[x],time_drift[y]))
  if index_d_yearday<>range(len(time_drift)):
 #     d_time=[d_time[i] for i in index_d_yearday]
      d_lat=[d_lat[i] for i in index_d_yearday]
      d_lon=[d_lon[i] for i in index_d_yearday]
      dep=[dep[i] for i in index_d_yearday]
      time_drift=sorted(time_drift)
      
  case=int(float(time_drift[-1]-time_drift[0])/numdays)#calculate case for each drifter
  

  
  startdate=num2date(time_drift[0])#give the startdate, in this example, startdate is start time of drift
  depth0=dep[0]
  print startdate  
  end_date=num2date(time_drift[-1])
  #Calculate the part of the data index for url
  for n in range(len(list_year)):
      #get part of the data from the url
      if list_year[n]==startdate.year:
          #for each year, get the index of start and end.
          index_start_url=list_year_index[n]
      if list_year[n]==end_date.year:
          index_end_url=list_year_index[n+1]


 # use function get_data_model to get the valuses from url
  part_url="["+str(index_start_url)+":1:"+str(index_end_url)+"]"
  url2='http://www.smast.umassd.edu:8080/thredds/dodsC/fvcom/hindcasts/30yr_gom3?'+'h[0:1:48450],lat[0:1:48450],lon[0:1:48450],nv[0:1:2][0:1:90414],siglay[0:1:44][0:1:48450],Times'+part_url+',u'+part_url+"[0:1:44][0:1:90414],"+"v"+part_url+"[0:1:44][0:1:90414]"
  jdmat_m,lat_vel_1,lon_vel_1,u,v,lat_vel,lon_vel,h_vel,siglay=get_FVCOM(url2)
  #plot the figure 
  fig=plt.figure(int(f))#,figsize=(16, 12))
  ax = fig.add_subplot(111)
  fig2=plt.figure(int(f+str(1)))
  ax2= fig2.add_subplot(111)
  
  
  lat_drift_all,lon_drift_all,=[],[],
  # for every two days, plot the drift and model track
  for k in range(0, case):#track_num):
    
    lat_drift_part,lon_drift_part=[],[]
    dist_drift,time_drift_id=[],[]
    #find the start date from drift valuses for each two days  
    index_start=int(round(np.interp(date2num(startdate)+2*k,time_drift,range(len(time_drift)))))
    
    if k<=1:
      print "index_start", index_start
    la=[d_lat[index_start]]
    lo=[d_lon[index_start]]

    startdate1=[num2date(time_drift[index_start])]
    depth0=dep[index_start]
    if k<=1:
      print "index_start", index_start,startdate1,
      print "timedrift", time_drift[0], startdate
      
    #use function gettrack to calculate the tracks  
    (lat_m,lon_m,time,uu,vv)=gettrack_FVCOM(depth0,jdmat_m,lat_vel_1,lon_vel_1,u,v,lat_vel,
                                  lon_vel,h_vel,siglay,startdate1[0],numdays,daystep,la[0],lo[0])#get the lat_model ,lon_model

    #use the startdate and starttime to calculate the drift's lat, lon 
    index_time=[]
    for i in time:
        index_time.append(np.interp(i,time_drift,range(len(time_drift))))
    
    for i in index_time:
        lat_drift_part.append(np.interp(i,range(len(d_lat)),d_lat))
        lon_drift_part.append(np.interp(i,range(len(d_lon)),d_lon))
        if k<3:
          lon_drift_all.append(np.interp(i,range(len(d_lon)),d_lon))
          lat_drift_all.append(np.interp(i,range(len(d_lat)),d_lat))

    if k<3:    
        ax.plot(lon_m,lat_m,"-",linewidth=1,marker=".")

    dist1=[]
    #calculate the separtcan distance
    for m in range(len(lon_m)):
        dist1.append(dist(lat_m[m],lon_m[m],lat_drift_part[m],lon_drift_part[m])[0])
        dist_drift.append(dist(lat_m[m],lon_m[m],lat_drift_part[m],lon_drift_part[m])[0])
        
    #calculate the each 2 day's mean distance
    two_days_mean_distance=np.mean(dist1)
    two_days_std_distance=np.std(dist1)
    if k <10:
        ax2.plot(time,dist1)
    locator = mpl.dates.AutoDateLocator()
    ax2.xaxis.set_major_locator(locator)
    monthsFmt = DateFormatter('%m/%d')
    ax2.xaxis.set_major_formatter(monthsFmt)
    ax2.set_xlabel(str(startdate1[0].year))
    ax2.set_ylabel("distance")
    
    if k<3:
      text1=ax.annotate(str(num2date(time[0]).month)+"-"+str(num2date(time[0]).day)+" " +str(num2date(time[0]).hour)+":"+str(num2date(time[0]).minute),
                      xy=(lon_m[0], lat_m[0]),  xycoords='data',size=15,
                      xytext=(5, 8), textcoords='offset points',
                      arrowprops=dict(arrowstyle="->"))

    time_0=time[0]
    for nn in time:
      time_drift_id.append(nn-time_0)
    ax3.plot(time_drift_id,dist_drift,'b')
    #for each certain hour, calculate the distance
    list_hour=list(scipy.linspace(0,numdays,48))
  

  
    for i in range(len(list_hour)):
        each_hour_days_dis.append(np.interp(list_hour[i],time_drift_id,dist_drift))
         
  ax.plot(lon_drift_all,lat_drift_all,"-",linewidth=3,marker=".",color="red")
  each_drift_std_distance=np.std(dist_drift)
  each_drift_mean_distance=np.mean(dist_drift)
 # fig.savefig(str(f)+"_driftVSmodel.png")
 # fig2.savefig(str(f)+"_timeVSds.png")
  fig.savefig("/net/nwebserver/epd/ocean/MainPage/fate/"+str(f)+"_driftVSmodel.png")
  fig2.savefig("/net/nwebserver/epd/ocean/MainPage/fate/"+str(f)+"_timeVSds.png")
  
  f_open.write('<tr><b><td align="right">' +'<a href="'+f+'_driftVSmodel.png'+'">'+f+'</a>'+'</td></b><b><td align="right">'+'<a href="'+str(f)+'_timeVSds'
               +'.png">'+' %7.3f'%each_drift_mean_distance+'</a>'+'</td></b>'
               '<b><td align="right"> %7.3f'%each_drift_std_distance+'</td></b>'+'<b><td align="right"> %7.0f'%case+'</td></b>'+'<b><td align="right"> %7.0f'%depth0+'</td></b>'+
               '</tr>\n')

  all_drift_std.append(each_drift_std_distance)
  all_drift_mean.append(each_drift_mean_distance)
#convert list to array and reshape it
each_hour_days_dis_array=scipy.array(each_hour_days_dis)
each_hour_dist=np.resize(each_hour_days_dis_array,(len(driftid),len(list_hour)))
dist_mean_all=[]
for k in range(0,len(list_hour)):
    each_hour_dist_all=[]
    for i in range(len(driftid)):
        each_hour_dist_all.append(each_hour_dist[i][k])
    dist_mean_all.append(np.mean(each_hour_dist_all))

ax3.plot(list_hour,dist_mean_all,"r",linewidth=3)
total_mean=np.mean(all_drift_mean)
total_std=np.mean(all_drift_std)
fig3.savefig("/net/nwebserver/epd/ocean/MainPage/fate/"+"sd_time.png")
#fig3.savefig("sd_time.png")
f_open.write('<tr><b><td align="right">' +"total"+'</td></b><b><td align="right">' +'<a href="sd_time.png">'+' %7.3f'%total_mean+'</a>'+'</td></b>'
               '<b><td align="right"> %7.3f'%total_std+'</td></b>'
               '</tr>\n')

f_open.close()
