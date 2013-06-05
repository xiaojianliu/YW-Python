from mpl_toolkits.mplot3d import Axes3D
import matplotlib.pyplot as plt
import numpy as np

fig = plt.figure()
ax = fig.add_subplot(111)

u = np.linspace(0, 2 * np.pi, 20)
v = np.linspace(0, np.pi,20)

x =  0.15*np.outer(np.cos(u), np.sin(v))
y =  0.15*np.outer(np.sin(u), np.sin(v))
#z =  np.outer(np.ones(np.size(u)), np.cos(v))
ax.plot(x[:,0:5]-67.28, y[:,0:5]+41.35,color='b')
plt.show()
X=list(x[:,0:5]-67.28)
Y=list(y[:,0:5]+41.35)
np.savez('XY',X=X,Y=Y)

#############above is 2D#######################

'''
fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')

u = np.linspace(0, 2 * np.pi,50)
v = np.linspace(0, np.pi,50)

x = np.outer(np.cos(u), np.sin(v))
y = np.outer(np.sin(u), np.sin(v))
z = np.outer(np.ones(np.size(u)), np.cos(v))
ax.plot_surface(x[:,0:25], y[:,0:25], z[:,0:25])

plt.show()
'''