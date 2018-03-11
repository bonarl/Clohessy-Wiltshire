#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Mar  9 16:22:46 2018

@author: bonar
"""

import math
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from matplotlib import animation

r0 = [5, 0.0, 0]                                #initial position vector
rdot0 = [0.1, 0.3, -0.023]                      #initial velocity
omeg = 0.0010854                                #mean motion of principal body, (can be calculated from altitude of circular orbit)

def CW(r0, rdot0, omeg, t):
    x0 = r0[0]
    y0 = r0[1]
    z0 = r0[2]
    xdot0 = rdot0[0]
    ydot0 = rdot0[1]
    zdot0 = rdot0[2]
    
    xt = (4*x0 + (2*ydot0)/omeg)+(xdot0/omeg)*math.sin(omeg*t)-(3*x0+(2*ydot0)/omeg)*math.cos(omeg*t)
    yt = (y0 - (2*xdot0)/omeg)+((2*xdot0)/omeg)*math.cos(omeg*t)+(6*x0 + (4*ydot0)/omeg)*math.sin(omeg*t)-(6*omeg*x0+3*ydot0)*t
    zt = z0*math.cos(omeg*t)+(zdot0/omeg)*math.sin(omeg*t)
    return([xt, yt, zt])
    
xs = []
ys = []
zs = []

nframes = 300
dt = 60
for j in range(nframes):
    xs.append(CW(r0, rdot0, omeg, j*dt)[0])
    ys.append(CW(r0, rdot0, omeg, j*dt)[1])
    zs.append(CW(r0, rdot0, omeg, j*dt)[2])
xmin = min(xs)
xmax = max(xs)
ymin = min(ys)
ymax = max(ys)
zmin = min(zs)
zmax = max(zs)

fig = plt.figure()
ax = fig.add_subplot(111, projection='3d',xlim = (xmin, xmax), ylim = (ymin, ymax),zlim = (zmin, zmax))
ax.set_xlabel('x')
ax.set_ylabel('y')
ax.set_zlabel('z')  
def init():
    return(ax)

def animate(j):
    x = CW(r0, rdot0, omeg, j*dt)[0]
    y = CW(r0, rdot0, omeg, j*dt)[1]
    z = CW(r0, rdot0, omeg, j*dt)[2]
    ax.scatter(x, y, z, marker='.',c='r', alpha = 0.5)
    return(ax)
    
anim= animation.FuncAnimation(fig, animate, init_func=init, frames = nframes, interval = 0.000001)
plt.show()
#anim.save('cw.gif', dpi=80, writer='imagemagick')
