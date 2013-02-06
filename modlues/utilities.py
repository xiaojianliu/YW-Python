
"""
Created on Thu Jan  5 14:40:30 2012
This is a collection of code I found outside normal Python distribution
that does a variety of things. For example, the "smooth" function below
does filtering better than I could find in the python moving average functions
@author: jmanning
"""

import numpy as np
import math
from datetime import timedelta as td
import matplotlib.dates as mdates
import colorsys

def choose_date_label(start,end):
          
    delta = end - start
    if delta <= td(minutes=10):
        loc = mdates.MinuteLocator()
        fmt = mdates.DateFormatter('%I:%M %p')
    elif delta <= td(minutes=30):
        loc = mdates.MinuteLocator(byminute=range(0,60,5))
        fmt = mdates.DateFormatter('%I:%M %p')
    elif delta <= td(hours=1):
        loc = mdates.MinuteLocator(byminute=range(0,60,15))
        fmt = mdates.DateFormatter('%I:%M %p')
    elif delta <= td(hours=6):
        loc = mdates.HourLocator()
        fmt = mdates.DateFormatter('%I:%M %p')
    elif delta <= td(days=1):
        loc = mdates.HourLocator(byhour=range(0,24,3))
        fmt = mdates.DateFormatter('%I:%M %p')
    elif delta <= td(days=3):
        loc = mdates.HourLocator(byhour=range(0,24,6))
        fmt = mdates.DateFormatter('%I:%M %p')
    elif delta <= td(weeks=2):
        loc = mdates.DayLocator()
        fmt = mdates.DateFormatter('%b %d')
    elif delta <= td(weeks=12):
        loc = mdates.WeekdayLocator()
        fmt = mdates.DateFormatter('%b %d')
    elif delta <= td(weeks=52):
        loc = mdates.MonthLocator()
        fmt = mdates.DateFormatter('%b')
    else:
        loc = mdates.MonthLocator(interval=3)
        fmt = mdates.DateFormatter('%b %Y')
    return loc,fmt

def fixticks(axis_object):
    print 'New2'  
    a = axis_object.get_xlim()
    loc,fmt = choose_date_label(a[0],a[1])
    axis_object.xaxis.set_minor_formatter(fmt)
    axis_object.xaxis.set_minor_locator(loc)


def lat2str(deg):
    min = 60 * (deg - np.floor(deg))
    deg = np.floor(deg)
    dir = 'N'
    if deg < 0:
        if min != 0.0:
            deg += 1.0
            min -= 60.0
        dir = 'S'
    return (u"%d\N{DEGREE SIGN} %g' %s") % (np.abs(deg),np.abs(min),dir)
  
def lon2str(deg):
    min = 60 * (deg - np.floor(deg))
    deg = np.floor(deg)
    dir = 'E'
    if deg < 0:
        if min != 0.0:
            deg += 1.0
            min -= 60.0
        dir = 'W'
    return (u"%d\N{DEGREE SIGN} %g' %s") % (np.abs(deg),np.abs(min),dir)

def nearxy(x,y,x0,y0):
    distance = []
    for i in range(0,len(x)):   
        distance.append(abs(math.sqrt((x[i]-x0)**2+(y[i]-y0)**2)))
    min_dis = min(distance)
    len_dis = len(distance)
    index = []
    for ii in range(0,len_dis):
        if distance[ii] == min_dis:
            index.append(ii)
    return min(distance),index[0]

def smooth(x,window_len=11,window='hanning'):
    """smooth the data using a window with requested size.
    
    This method is based on the convolution of a scaled window with the signal.
    The signal is prepared by introducing reflected copies of the signal 
    (with the window size) in both ends so that transient parts are minimized
    in the begining and end part of the output signal.
    
    input:
        x: the input signal 
        window_len: the dimension of the smoothing window; should be an odd integer
        window: the type of window from 'flat', 'hanning', 'hamming', 'bartlett', 'blackman'
            flat window will produce a moving average smoothing.

    output:
        the smoothed signal
        
    example:

    t=linspace(-2,2,0.1)
    x=sin(t)+randn(len(t))*0.1
    y=smooth(x)
    
    see also: 
    
    numpy.hanning, numpy.hamming, numpy.bartlett, numpy.blackman, numpy.convolve
    scipy.signal.lfilter
 
    TODO: the window parameter could be the window itself if an array instead of a string   
    """

    if x.ndim != 1:
        raise ValueError, "smooth only accepts 1 dimension arrays."

    if x.size < window_len:
        raise ValueError, "Input vector needs to be bigger than window size."

    if window_len < 3:
        return x

    if not window in ['flat', 'hanning', 'hamming', 'bartlett', 'blackman']:
        raise ValueError, "Window is on of 'flat', 'hanning', 'hamming', 'bartlett', 'blackman'"

    s = np.r_[x[window_len-1:0:-1], x, x[-1:-window_len:-1]]
    
    if window == 'flat': #moving average
        w = np.ones(window_len,'d')
    else:
        w = eval('numpy.'+ window +'(window_len)')

    y = np.convolve(w/w.sum(),s,mode='valid')
    return y


def uniquecolors(N):
    """
    Generate unique RGB colors
    input: number of RGB tuples to generate
    output: list of RGB tuples
    """
    HSV_tuples = [(x*1.0/N, 0.5, 0.5) for x in range(N)]
    RGB_tuples = map(lambda x: colorsys.hsv_to_rgb(*x), HSV_tuples)
    colors =  [colorsys.hsv_to_rgb(x*1.0/N, 0.5, 0.5) for x in range(N)]
    return colors

def colors(n):
	  """Compute a list of distinct colors, each of which is represented as an RGB 3-tuple."""
	  """It's useful for less than 100 numbers"""
	  if pow(n,float(1)/3)%1==0.0:
	     n+=1 
	  #make sure number we get is more than we need.
	  rgbcolors=[]
	  x=pow(n,float(1)/3)
	  a=int(x)
	  b=int(x)
	  c=int(x)
	  if a*b*c<=n:
	    a+=1
	  if a*b*c<n:
	    b+=1
	  if a*b*c<n:
	    c+=1
	  for i in range(a):
	      r=0.99/(a)*(i)
	      for j in range(b):
	          s=0.99/(b)*(j)
	          for k in range(c):
	              t=0.99/(c)*(k)
	              color=r,s,t
	              rgbcolors.append(color)
	  return rgbcolors
   


   