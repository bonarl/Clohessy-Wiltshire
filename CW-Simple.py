### https://github.com/bonarl/Clohessy-Wiltshire/blob/master/cw.py

import math
import matplotlib.pyplot as plt
from matplotlib import animation
from mpl_toolkits.mplot3d import Axes3D

r0    = [5.000, 0.000,  0.000]   # initial position vector
rdot0 = [0.100, 0.300, -0.023]   # initial velocity
omega =  0.0010854               # mean motion of principal body (calculated from circular orbit altitude)

def CW(r0, rdot0, omega, t):
  x0    = r0[0]
  y0    = r0[1]
  z0    = r0[2]
  xdot0 = rdot0[0]
  ydot0 = rdot0[1]
  zdot0 = rdot0[2]

  A = [[         4.0 - 3.0*math.cos(omega*t), 0.0,               0.0,           math.sin(omega*t)/omega,   2.0*(1.0-math.cos(omega*t))/omega,                     0.0],
       [ 6.0*math.sin(omega*t) - 6.0*omega*t, 1.0,               0.0, 2.0*(math.cos(omega*t)-1.0)/omega, 4.0*math.sin(omega*t)/omega - 3.0*t,                     0.0],
       [                                 0.0, 0.0, math.cos(omega*t),                               0.0,                                 0.0, math.sin(omega*t)/omega]]
  
  xt = A[0][0]*r0[0] + A[0][1]*r0[1] + A[0][2]*r0[2] + A[0][3]*rdot0[0] + A[0][4]*rdot0[1] + A[0][5]*rdot0[2]
  yt = A[1][0]*r0[0] + A[1][1]*r0[1] + A[1][2]*r0[2] + A[1][3]*rdot0[0] + A[1][4]*rdot0[1] + A[1][5]*rdot0[2]
  zt = A[2][0]*r0[0] + A[2][1]*r0[1] + A[2][2]*r0[2] + A[2][3]*rdot0[0] + A[2][4]*rdot0[1] + A[2][5]*rdot0[2]
  return([xt, yt, zt])

xs = []
ys = []
zs = []

nframes = 300
dt = 60

for j in range(nframes):
  xs.append(CW(r0, rdot0, omega, j*dt)[0])
  ys.append(CW(r0, rdot0, omega, j*dt)[1])
  zs.append(CW(r0, rdot0, omega, j*dt)[2])
xmin = min(xs)
xmax = max(xs)
ymin = min(ys)
ymax = max(ys)
zmin = min(zs)
zmax = max(zs)

fig = plt.figure()
ax = fig.add_subplot(111, projection = '3d', xlim = (xmin, xmax), ylim = (ymin, ymax), zlim = (zmin, zmax))
ax.set_xlabel('x')
ax.set_ylabel('y')
ax.set_zlabel('z')  

def init():
  return(ax)

def animate(j):
  x = CW(r0, rdot0, omega, j*dt)[0]
  y = CW(r0, rdot0, omega, j*dt)[1]
  z = CW(r0, rdot0, omega, j*dt)[2]
  ax.scatter(x, y, z, marker = '.', c = 'r', alpha = 0.5)
  return(ax)

anim = animation.FuncAnimation(fig, animate, init_func=init, frames = nframes, interval = 0.001)
anim.save('cw.gif', dpi=80, writer='Pillow')