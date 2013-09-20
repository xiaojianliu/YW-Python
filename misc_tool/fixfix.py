################################################################################################
## fixfix: post-processes drifter tracks
## xiuling.wu
## August 9, 2011
##
## Jim Manning modifications
## October 2011
## March 2012
## Nov 2012  trouble with ginput causing segfault
## STEPS ARE AS FOLLOWS:
#  - load all the ascii data file
#  - for each drifter it conducts 4 basic steps
#    - eliminate repeat times
#    - calculate forward and backward differences (velocity) and eliminate bad points
#    - writes meta data to log file
#    - check for drogue loss
#    - clean plots of track (pth) and velocity uv_id_final.png)
#  - generates oracle ready ascii file
##################################################################################################

import sys
sys.path.append("/home3/ocn/jmanning/py/jmanning/whython6/")
import csv
from  conversions import ll2uv #To convert from longitude/latitude to unit vectors
import numpy as np
import matplotlib as mpl
import matplotlib.mlab as ml
#part1
#import scipy
#from scipy import *
#import pylab
from pylab import *
import matplotlib.pyplot as plt
#import basemap
from matplotlib.dates import num2date,date2num, DateFormatter
import math


### HARD CODE ################################################################################

critfactor=8 #multipe of the mean velcity to discard
minla=20;maxla=48;minlo=-150;maxlo=-20  # delimiting data to make it easier
bathy=True # set to "True" if isobaths are wanted

fid=open("/home3/ocn/jmanning/drift/drift.log","a")#permanent log file that is appended to../home3/ocn/jmanning/drift/drift.log

# load /data5/jmanning/drift/massdmf/withtemp_2009.dat   ##ids with temperature
#list(set(line.strip() for line in open('X:/drift/bowdoin/2009/withtemp_2009.dat'))) # ids with temperature
#fileformat='bowdoin'; # format of temperature data (sometimes "emolt") if using minilog/dat files
#wt='withtemp_2009'

direcin="/net/nwebserver/drifter/" # directory where the final plots are stored
#direcout='/net/data5/jmanning/drift/umassb/2012'  #'/data5/jmanning/drift/css/2011'
#fn='drift_umb_2012_1.dat'
direcout='/net/data5/jmanning/drift/kocik/2012'  #'/data5/jmanning/drift/css/2011'
fn='drift_cfcc_2012_1.dat'
fid.write('\n'+'#'*40+' below is '+str(fn)+' log '+'#'*40+'\n')
depcont=[-10] # depth contour in final plot (apparently not used in 5/2012)
year=int(2012)
strattle=0 # this will eventually flag to "1" if yeardays change by more than 365
### END HARDCODES ################################################

#raw data loaded
idss,yeardays_all,lat_all,lon_all,day_all,hours_all,minutes_all,depths,temps=[],[],[],[],[],[],[],[],[]
csv_reader=csv.reader(open(direcin+fn,"r"))

for line in (x for x in csv_reader if x[0][0] !='%'):  # if the first line is comment line, skip
  #  print float(line[0].split()[8])
    if float(line[0].split()[8])<maxla and float(line[0].split()[8])>minla and float(line[0].split()[7])>minlo and float(line[0].split()[7])<maxlo:
        idss.append(line[0].split()[0])
        yeardays_all.append(line[0].split()[6])
        lat_all.append(line[0].split()[8])
        lon_all.append(line[0].split()[7])
        day_all.append(line[0].split()[3])
        hours_all.append(line[0].split()[4])
        minutes_all.append(line[0].split()[5])
        depths.append(line[0].split()[9])
        temps.append(line[0].split()[10])

# get ids
id=list(set(idss))

# convert string to float
yeardays_all=[float(i)+1 for i in yeardays_all]# in python num2date(), less one day than matlab, so add 1 here
lat_all=[float(i) for i in lat_all]
lon_all=[float(i) for i in lon_all]

#days_all=[float(i) for i in days_all]
#ids=[float(i) for i in ids]
days_all,ids=[],[]
for i in range(len(day_all)):
    days_all.append(int(float(day_all[i])))
for i in range(len(id)):
    ids.append(int(float(id[i])))

  
fido=open(direcout+'/prep_for_oracle_'+fn[6:],'w')

ids=np.sort(ids)
#ids=[110410712]#,11041073]
for k in range(len(ids)): #where "ids" is a list of distinct ids and int
    # latitude, longitude, time
    strattle=0
    lat,lon,time,yeardays,depth,temp=[],[],[],[],[],[]
    for i in range(len(idss)):
        if int(float(idss[i]))==ids[k]:
            lat.append(lat_all[i])
            lon.append(lon_all[i])
            if (i>1)&(strattle==0):# here's where we account for a new year
                if yeardays_all[i]-yeardays_all[i-1]<-200:
                   year=year+1
                   print 'incremented year to '+str(year)
                   strattle=1
            yeardays.append(yeardays_all[i])
            #time.append(date2num(num2date(yeardays_all[i]).replace(year=year).replace(day=days_all[i])))
            time.append(date2num(num2date(yeardays_all[i]).replace(year=year)))
            depth.append(depths[i])
            temp.append(temps[i])
    #print time        
    print "there are ", len(lat), " fixes for id =",ids[k]          
    print "Note: Check to see if any of these already exist in database before loading"

    # STEP 1a: check to make sure time is monotonic
    #if len(find(diff(time)<=0))>0:
    #  plot(time)
    #  show()
      #raw_input('Trouble with time not increasing press return to continue')
    #  close()
      
    # STEP 1b: check for repeat times
    ###### while time[i]==time[i-1], get the del_same_time_index ########
    del_same_time_index=[]
    for i in range(1,len(time)):
        if int(time[i-1])==int(time[i]) and num2date(time[i-1]).hour== num2date(time[i]).hour and  num2date(time[i-1]).minute== num2date(time[i]).minute:
            del_same_time_index.append(i)
    del_same_time_index.reverse()
    if del_same_time_index==[]:
        print "there is no same time."
    
    else:
        print str(len(del_same_time_index))+' points deleted with the same times'
        index=range(len(time))
        for i in del_same_time_index:
            del lat[i],lon[i],time[i],yeardays[i],depth[i],temp[i]
            
    # STEP 2a:
    ############ calculate forward and backward velocities of the raw data  ##########################

        
    forward_u,forward_v,forward_spd,jdn=ll2uv(time,lat,lon)# was yeardays but now uses "time" as of 3/29/2012
    backward_u,backward_v,backward_spd,backward_jdn=ll2uv(time[::-1],lat[::-1],lon[::-1])

    ## calculate resultants
    id_fs=list(np.where(np.array(forward_spd)<500)[0])
    id_bs=list(np.where(np.array(backward_spd)<500)[0])
    idc=[val for val in id_fs if val in id_bs]
    jdraw,spdraw=[],[]
    for i in idc:
        jdraw.append(jdn[i])
        spdraw.append(forward_spd[i])
    ###########  plot the velocities  ###################################################
##    def plot_speed(time,speed):
##        #fig=plt.figure()
##        ax = fig.add_subplot(111) #to divide the fig into some area and (line row chose) 
##        plt.title('Difter#'+str(ids[k]))
##        jd=[]
##        for i in time:
##            jd.append(i)
##        plt.plot(jd,speed,'b-')
        #locator = mpl.dates.AutoDateLocator()
        #ax.xaxis.set_major_locator(locator)
##        if len(jd)<100:
##        else:
##            monthsFmt = DateFormatter('%b/%d')
##        ax.set_ylabel('cm/s')
        #ax.xaxis.set_major_formatter(monthsFmt)
##        ax.set_xlabel(str(year),fontsize=17)
##        plt.grid()
##    fig=plt.figure()
##    plot_speed(jdraw,spdraw)## plot speed
##    plt.show()
##    plt.close()

    #######################################################################################
    
    # calculate a reasonable criteria for this drifter
    crit=np.mean([abs(i) for i in forward_spd])*critfactor
    print "Velocity criteria set to ", str(critfactor),' times the mean or ',str(crit),' cm/s'
    # check for drifter going aground (ie very low velocity)
    idlow=list(np.where(np.array(spdraw)<float(np.mean([abs(i) for i in forward_spd]))/100)[0])
    # if idlow is not empty, add the comments in fid file
    if idlow<>[]:
        for i in range(len(idlow)):
            print 'WARNING: Drifter ',str(ids[k]),' may be hung up on gear or aground on ',str(idlow[i]),' where velocity is < 1# mean'
            #fid.write(str(ids[k]).rjust(10)+' apparently hung-up on '+str(idlow[i])+'\n')
        idlow_print0=str(sorted(idlow))
        idlow_print1=idlow_print0.replace(', ',' ')
        tempochunk0=str(ids[k]).rjust(10)+' apparently hung-up on '+str(idlow_print1)+'\n'#'from'+str(idlow[0])+'to'+str(idlow[-1])+'\n'
    else:
        tempochunk0='There is no point hung up\n'


    #### find bad velocities where criteria was just calculated
    idbadf=list(np.where(abs(np.array(forward_spd))>crit)[0])
    idbadb=list(np.where(abs(np.array(backward_spd))>crit)[0])
    
    #if it is the 2nd time/point in the bad forward velocity (fu) that caused the problem
    # then the 2nd time/point associated with the bad backward velocity should match
    timeb=time[::-1] # backwards time vector
    badtime=list(set([time[i+1] for i in idbadf]).intersection(set([timeb[i+1] for i in idbadb])))
    print "%10.3f percent bad velocities deleted according to velocity criteria" % float(len(badtime)/float(len(lat))*100.)
    index_badtime=[]# find near the badtime points
    for i in badtime:
        index_badtime.append(int(np.interp(i,time,range(len(time)))))
    if index_badtime<>[]:
      index_badtime.reverse()
      for i in index_badtime:
        index_near_badtimes=[]
        if i-5<0: 
            ra=range(0,i+5)
        elif i+5>len(lat):
            ra=range(i-5,len(lat)-1)
        else:
            ra=range(i-5,i+5)
        for m in ra:
            index_near_badtimes.append(m)
        plot_badtime=list(set(index_near_badtimes))
        #plot the bad time data and near the bad time data's points
        #plt.plot([lon[l] for l in plot_badtime],[lat[l] for l in plot_badtime],marker='.',)
        #plt.plot([lon[l] for l in index_badtime],[lat[l] for l in index_badtime],marker='o',markerfacecolor='r',linestyle='None')
        fig=plt.figure() 
        plt.plot([lon[l] for l in plot_badtime],[lat[l] for l in plot_badtime],marker='.',)
        plt.plot(lon[i],lat[i],marker='o',markerfacecolor='r',linestyle='None')
        plt.show()
        #plt.close()
        del_or_not=raw_input('Delete? (y/n or 1 for the end point): ')
        if del_or_not=='y':
           del time[i],lat[i],lon[i],yeardays[i],depth[i],temp[i]
        elif del_or_not=='1':
           plt.plot(lon[i-1],lat[i-1],marker='o',markerfacecolor='r',linestyle='None')
           plt.show()
           raw_input('How is that? press return')
           del time[i-1],lat[i-1],lon[i-1],yeardays[i-1],depth[i-1],temp[i-1]
           plt.close()
        plt.close()  
    
    # STEP 3:
    # delelte values bad due to objective criteria    
    #index_badtime.reverse()
    #for i in index_badtime:
    #    if lat[i]!=999:
    #      del time[i],lat[i],lon[i],yeardays[i],depth[i],temp[i]
    #print str(float(len(badtime))/len(time)*100),'# editted due to bad velocities > criteria'
    idgood=len(lat)
    ###############################################################################################

    
    # Step 4a:
    # calculate forward velocities of the automatically editted data

    fu,fv,spd1,jd1=ll2uv(time,lat,lon)
    fig=plt.figure()
    #plot_speed(jd1,spd1)
    plt.plot(jd1,spd1)
    plt.plot(jd1,spd1,marker="o",markerfacecolor="r",linestyle='None')
    plt.show()
    print 'click on any obviously bad points and then press the enter key.'
    badpoints=ginput(n=0)
    print badpoints#,timeout=10)
    #badpoints=ginput(0,timeout=10,mouse_stop=3)
   
    #badpoints=[]
    #close()
    
    # Step 4b:
    # eliminate those points clicked as bads

    # find badpoints index in yeardays
    index_badpoints=[]
    badpoints_num=len(badpoints)
    for i in range(len(badpoints)):
        index_badpoints.append(int(np.interp(badpoints[i][0],jd1,range(len(jd1)))))
    print index_badpoints
    index_badpoints=list(set(index_badpoints))
    print "%10.2f percent bad velocities deleted according to manual clicks on velocity" % float(float(badpoints_num)/len(lat)*100.)
    for i in sorted(index_badpoints)[::-1]:
        del time[i], lat[i], lon[i], yeardays[i],depth[i],temp[i]
    
    #plt.close()
   
      
        
    plot_again=raw_input("Do you want to replot the figure after delete the bad points?(y/n)")
    #plot_again='y'
    if plot_again=="y" or plot_again=="Y" or plot_again=="yes":
        #plt.close('all')
        fig=plt.figure()
        fu2,fv2,spd2,jd2=ll2uv(time,lat,lon)
        #plot_speed(jd2,spd2)
        plt.plot(jd2,spd2,'mo-',markersize=5)#marker='o',markerfacecolor="r",linestyle='None')
        plt.show()
        #plt.close()
    #print 'pausing 10 seconds'
    #sleep(10)
    
    # if there are a list of bad points, click the first point and the last point, then delete between them
    del_between=raw_input('Do you want to delete all the points between two points? input "N" or "Y"' )
    #del_between='N'
    if del_between=="N" or del_between=="n":
        print "You have choosen NOT to delete all the points between two points."
    if del_between=="Y" or del_between=="y" :
        print "Please click the first bad point and the last bad point to choose the range of the bad points"
        between_badpoints=ginput(n=0)
        print between_badpoints#,timeout=0)#,mouse_stop=2)
        index_between_badpoints=[]
        for i in range(len(between_badpoints)):
            index_between_badpoints.append(int(np.interp(between_badpoints[i][0],jd2,range(len(jd2)))))
        print index_between_badpoints
        index_betweens=[]
        for i in range(sorted(index_between_badpoints)[0],sorted(index_between_badpoints)[1]+1):
            index_betweens.append(i)
        for i in index_betweens[::-1]:
            del lat[i],lon[i],time[i],yeardays[i],depth[i],temp[i]
        del_between_badpoints=sorted(index_between_badpoints)[1]-sorted(index_between_badpoints)[0]+1
        badpoints_num=len(badpoints)+del_between_badpoints
        print "%10.2f percent editted due to bad velocities from manual clicks between two points" % float(float(badpoints_num)/len(time)*100.)
    
    if ids[k]==1174306915 or ids[k]==1174306911:
        del time[-1], lat[-1], lon[-1], yeardays[-1],depth[-1],temp[-1] 
        
    fig=plt.figure()
    fu3,fv3,spd3,jd3=ll2uv(time,lat,lon)
    #plot_speed(jd3,spd3)
    plt.plot(jd3,spd3,'bo-')  
    plt.show()

    #step 5a:
    #manually delete points based on track
 ##############################################################################################################   
        
    fig=plt.figure()
    #plt.figure(2)
    #plt.plot(lon,lat,'k-')
    plt.plot(lon,lat,'ro-')
    plt.show()
    print 'click on any obviously bad points and then press the enter key on the track.'
    bad=ginput(n=0)
    print bad
    badplotpts=[] #index of points that are found to be near the same index of x & y
    if len(bad)>0:
      for kbad in range(len(bad)):
        idxbad=ml.find(abs(lon-bad[kbad][0])==min(abs(lon-bad[kbad][0])))
        idybad=ml.find(abs(lat-bad[kbad][1])==min(abs(lat-bad[kbad][1])))
        print idxbad,idybad
        if idxbad==idybad:
          print lat[int(idxbad)],lon[int(idxbad)],' is bad'  
          badplotpts.append(int(idxbad))    
      for kk in range(len(badplotpts)):
        plt.plot(lon[badplotpts[kk]],lat[badplotpts[kk]],'bo')
#      #thismanager = get_current_fig_manager()
#      #thismanager.window.SetPosition((1000,500))
        plt.show()
      for i in sorted(badplotpts)[::-1]:
        del time[i], lat[i], lon[i], yeardays[i],depth[i],temp[i]
      
      #fig=plt.figure()
      plt.plot(lon,lat,'yo-')
      plt.show()
      raw_input(str(len(badplotpts))+' deleted from manual click on track. Press return to continue')  
        #plt.close()
   
    
    # write to log file if some data was editted
    if badpoints_num>0:
        tempochunk1=(str(ids[k]).rjust(10)+' '+ str(crit).rjust(10)+' '+ str(badpoints_num).rjust(10)+' '+
                  str(idgood).rjust(10)+" "+str(math.floor(time[-1]-time[0])).rjust(10)+" manual editted uv plot\n")
    else:
        tempochunk1='There is no bad point delete manual.\n'
    if len(badtime)>0:
        tempochunk2=(str(ids[k]).rjust(10)+' '+ str(crit).rjust(10)+' '+ str(len(badtime)).rjust(10)+' '+
                  str(idgood).rjust(10)+" "+str(math.floor(time[-1]-time[0])).rjust(10)+" objectively editted\n")
    else:
        tempochunk2='There is no bad velocities deleted according to velocity criteria.\n'
              
    if len(badplotpts)>0:
        tempochunk3=(str(ids[k]).rjust(10)+' '+ str(crit).rjust(10)+' '+ str(len(badplotpts)).rjust(10)+' '+
                  str(idgood).rjust(10)+" "+str(math.floor(time[-1]-time[0])).rjust(10)+" manually editted track points\n")
    else:
        tempochunk3='There is no bad point delete manual on track.\n'

    # clean velocities w/out bad points
    [u2,v2,spd2,jd2]=ll2uv(time,lat,lon) #was "yeardays" until 3/2012 to deal with strattling new year

    #plot time, lat,lon
    fig=plt.figure(figsize=(10,8))
    ax = fig.add_subplot(111)
    print 'calling basemap_usgs '
    #print min(lat),min(lon),max(lat),max(lon)
    #if max(lon)-min(lon)>4.0:
    #  basemap.basemap(lat,lon)#,(float(min(max(lat)-min(lat),max(lon)-min(lon)))+1.0)/5*4)
    #else:  
    #  basemap.basemap_detail(lat,lon,bathy, False,float(min(max(lat)-min(lat),max(lon)-min(lon)))/5*4)
    #basemap.basemap_usgs(lat,lon,False)#,depcont)
    print 'past basemap'
    ax.plot(lon,lat,marker=".",markerfacecolor='r',markersize=10)
    points_num=10
    ax.set_ylabel('latitude')
    ax.set_xlabel("longitude")
    #ax.set_xlim(min(lon)-(max(lon)-min(lon))/10.,max(lon)+(max(lon)-min(lon))/10.)
    #ax.set_ylim(min(lat)-(max(lat)-min(lat))/10.,max(lat)+(max(lat)-min(lat))/10.)
    ax.set_xlim(min(lon),max(lon))
    ax.set_ylim(min(lat),max(lat))
    plt.title('Drifter #'+str(ids[k]),fontsize=16)
    #  fig.autofmt_xdate()  #display the time "lean"
    ax.xaxis.set_label_coords(0.5, -0.1)#set the position of the xtick labels
    ax.yaxis.set_label_coords(-0.08, 0.4)
    # the last points, annotate "end"
    self_annotate1=ax.annotate("End", xy=(lon[-1], lat[-1]),xycoords='data', xytext=(8, 11), 
                              textcoords='offset points',arrowprops=dict(arrowstyle="->"))
    self_annotate2=ax.annotate("Start", xy=(lon[0], lat[0]),xycoords='data', xytext=(8, 11), 
                              textcoords='offset points',arrowprops=dict(arrowstyle="->"))
    if time[-1]-time[0]<=2:
        if len(time)<5: skip=1
        else: skip=int(float(len(time))/5)
        for i in range(0,len(time),skip):
            self_annotate3=ax.annotate(num2date(time[i]).replace(tzinfo=None).strftime('%d-%b %H:%M'), xy=(lon[i], lat[i]),
                                      xycoords='data', xytext=(8, 11), textcoords='offset points',arrowprops=dict(arrowstyle="->"))
       #     self_annotate3.draggable()
    elif (time[-1]-time[0]>2.0)&(time[-1]-time[0]<20.0):
        for i in range(1,len(time)):
            if num2date(time[i-1]).day<>num2date(time[i]).day:
                self_annotate4=ax.annotate(num2date(time[i]).replace(tzinfo=None).strftime('%d-%b'), xy=(lon[i], lat[i]),
                                          xycoords='data', xytext=(8, 11), textcoords='offset points',arrowprops=dict(arrowstyle="->"))
    else: # place approximately 10 labels
        for i in range(1,len(time),int(len(time)/10.)):
           #if num2date(time[i-1]).day<>num2date(time[i]).day:
           self_annotate4=ax.annotate(num2date(time[i]).replace(tzinfo=None).strftime('%d-%b'), xy=(lon[i], lat[i]),
                                          xycoords='data', xytext=(8, 11), textcoords='offset points',arrowprops=dict(arrowstyle="->"))          #      self_annotate4.draggable() #drag the text if you want
       

    thismanager = plt.get_current_fig_manager()
    thismanager.window.SetPosition((1000, 0))
    plt.show()
    plt.savefig(direcin+'pth_'+str(ids[k])+'_final'+".ps")
    plt.savefig(direcin+'pth_'+str(ids[k])+'_final'+".png")
    raw_input('press return to close final track window')
    plt.close()
   # plt.show()

    #plot u & v
    fig=plt.figure()
    ax1 = fig.add_subplot(111)
    plt.plot(jdn,forward_u,"r",label='raw eastward')
    plt.plot(jdn, forward_v,"b",label='raw northward')
    plt.plot(jd2,u2,"m",linewidth=2,label='final eastward')
    plt.plot(jd2,v2,"g",linewidth=2,label='final northward')
    leg=plt.legend()
   # leg.draggable()

    locator = mpl.dates.AutoDateLocator()
    ax1.xaxis.set_major_locator(locator)
    if len(jdn)<100:
        monthsFmt = DateFormatter('%b/%d %H:')
    else:
        monthsFmt = DateFormatter('%b/%d')
    ax1.xaxis.set_major_formatter(monthsFmt)
    ax1.set_xlabel(str(year))
    ax1.set_ylabel('cm/s (where 50 ~ 1 knot)')
    fig.autofmt_xdate()  #display the time "lean"
    plt.title('Drifter '+str(ids[k])+' cleaned',fontsize=16)
    plt.savefig(direcin+'uv_'+str(ids[k])+'_final'+'.ps')
    plt.savefig(direcin+'uv_'+str(ids[k])+'_final'+'.png')
    plt.show()
    raw_input('press return to close uv window')
   # close()

    # write out id,date,lat,lon,yrday0_gmt,temp, and depth_i

    depth=[float(i) for i in depth]
    for i in range(len(time)):     
        fido.write(str(ids[k]).rjust(10)+ " " +num2date(time[i]).replace(tzinfo=None).strftime('%d-%b-%Y:%H:%M')+" ")
        fido.write(("%10.6f") %(lat[i]))
        fido.write(" ")
        fido.write(("%10.6f") %(lon[i]))
        fido.write(" ")
        fido.write(("%10.6f") %(yeardays[i]-1))
        fido.write(" ")
        fido.write(temp[i]+ " ")
        fido.write(("%5.1f") %(depth[i]))
        fido.write('\n')
    if k<>len(ids)-1:        
      raw_input("Press Enter to process next drifter")  

    whetherlog=raw_input('Do you want to keep this log?')
    if whetherlog=="Y" or whetherlog=="y" :
        fid.write(tempochunk0)
        fid.write(tempochunk1)
        fid.write(tempochunk2)
        fid.write(tempochunk3)
    print 'log has been saved.'
fido.close()
fid.close()

