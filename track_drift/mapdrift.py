import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.basemap import Basemap
latfile='20010129latt.npy'
lonfile='20010129lont.npy'
lat=np.load(latfile)
lon=np.load(lonfile)
latsize=[min(lat)-2,max(lat)+2]
lonsize=[min(lon)-2,max(lon)+2]
plt.figure()
m = Basemap(projection='cyl',llcrnrlat=min(latsize),urcrnrlat=max(latsize),\
            llcrnrlon=min(lonsize),urcrnrlon=max(lonsize),resolution='h')#,fix_aspect=False)
m.drawparallels(np.arange(int(min(latsize)),int(max(latsize))+1,1),labels=[1,0,0,0])
m.drawmeridians(np.arange(int(min(lonsize)),int(max(lonsize))+1,1),labels=[0,0,0,1])
m.drawcoastlines()
m.fillcontinents(color='grey')
m.drawmapboundary()
m.plot(lon,lat,linewidth=1.5,color='r')
plt.title(latfile[:8]+' drift track')
plt.savefig(latfile[:8]+'Drift_track.png')
plt.show()