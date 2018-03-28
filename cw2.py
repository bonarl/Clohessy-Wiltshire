#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Mar 23 13:38:54 2018

@author: bonar
"""

#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Mar  9 16:22:46 2018

@author: bonar
"""

import math
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

plt.style.use('fivethirtyeight')

r0 = [0, 0, 0]                  #initial starting position
rdot0 = [0.1, 0, 0]            #initial velocity (on release from ISS) (m/s)
R = 405 + 6870   
mu = 398600.5
omeg = math.sqrt(mu/R**3)

nframes = 1000
nframes2 = 200
nframes3 = 200
dt = 20

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
def CW2(r0, rdot0, omeg, t):
    x0 = r0[0]
    y0 = r0[1]
    z0 = r0[2]
    xdot0 = rdot0[0]
    ydot0 = rdot0[1]
    zdot0 = rdot0[2]
    
    xt = (4*x0 + (2*ydot0)/omeg)+(xdot0/omeg)*math.sin(omeg*t)-(3*x0+(2*ydot0)/omeg)*math.cos(omeg*t)
    yt = (y0 - (2*xdot0)/omeg)+((2*xdot0)/omeg)*math.cos(omeg*t)+(6*x0 + (4*ydot0)/omeg)*math.sin(omeg*t)-(6*omeg*x0+3*ydot0)*t
    zt = z0*math.cos(omeg*t)+(zdot0/omeg)*math.sin(omeg*t)
    
    xdott = (3*omeg*x0+2*ydot0)*math.sin(omeg*t)+xdot0*math.cos(omeg*t)
    ydott = (6*omeg*x0+4*ydot0)*math.cos(omeg*t)-2*xdot0*math.sin(omeg*t)-(6*omeg*x0+3*ydot0)
    zdott = zdot0*math.cos(omeg*t)-z0*omeg*math.sin(omeg*t)
    
    return([xt,yt,zt],[xdott,ydott,zdott])
        
xs = []
ys = []
zs = []
ds = []



d =0
d_max = 50
stindex = nframes
for i in range(nframes):
    t = dt*i
    if d<d_max:
        r_vec= CW2(r0, rdot0, omeg, t)[0]
        rdot_vec = CW2(r0, rdot0, omeg, t)[1]
        x = r_vec[0]
        y = r_vec[1]
        z = r_vec[2]
        xdot = rdot_vec[0]
        ydot = rdot_vec[1]
        zdot = rdot_vec[2]
        rdot = [xdot, ydot, zdot]
        d = math.sqrt(x**2+y**2+z**2)
        ds.append(d)
        
        xs.append(x)
        ys.append(y)
        zs.append(z)
        
        v = math.sqrt(xdot**2+ydot**2+zdot**2)
        stindex = i
    
#second manouvre
     
r0 = [x, y, z]
deltav = [-0.03, 0, 0.02]
for i in range(len(rdot)):
    rdot0[i] = rdot[i]+deltav[i]


for i in range(nframes2):
    t = dt*i
    r_vec= CW2(r0, rdot0, omeg, t)[0]
    rdot_vec = CW2(r0, rdot0, omeg, t)[1]
    x = r_vec[0]
    y = r_vec[1]
    z = r_vec[2]
    xdot = rdot_vec[0]
    ydot = rdot_vec[1]
    zdot = rdot_vec[2]
    rdot = [xdot, ydot, zdot]
    d = math.sqrt(x**2+y**2+z**2)
    ds.append(d)
    
    xs.append(x)
    ys.append(y)
    zs.append(z)
    
    v = math.sqrt(xdot**2+ydot**2+zdot**2)

#maonoeuvre 3

r0 = [x, y, z]
deltav = [0, -0.05, 0.03]
for i in range(len(rdot)):
    rdot0[i] = rdot[i]+deltav[i]  
  

for i in range(nframes3):
    t = dt*i
    r_vec= CW2(r0, rdot0, omeg, t)[0]
    rdot_vec = CW2(r0, rdot0, omeg, t)[1]
    x = r_vec[0]
    y = r_vec[1]
    z = r_vec[2]
    xdot = rdot_vec[0]
    ydot = rdot_vec[1]
    zdot = rdot_vec[2]
    
    d = math.sqrt(x**2+y**2+z**2)
    ds.append(d)
    
    xs.append(x)
    ys.append(y)
    zs.append(z)
    
    v = math.sqrt(xdot**2+ydot**2+zdot**2)

xmin = min(xs)
xmax = max(xs)
ymin = min(ys)
ymax = max(ys)
zmin = min(zs)
zmax = max(zs)
rmax = max(ds)




fig = plt.figure()
ax = fig.add_subplot(111, projection='3d',xlim = (xmin, xmax), ylim = (ymin, ymax),zlim = (zmin, zmax))
ax.set_xlabel('x')
ax.set_ylabel('y')
ax.set_zlabel('z')  


frame = 0
while frame < stindex:
    frame +=1
    ax.scatter(xs[frame], ys[frame], zs[frame], marker = '.', color = 'g', alpha = 0.5, s=1)
    
while frame < nframes2+stindex:
    frame += 1
    ax.scatter(xs[frame], ys[frame], zs[frame], marker = '.', color = 'skyblue', s=1)
    
while frame < nframes3+nframes2+stindex:
    frame += 1
    ax.scatter(xs[frame], ys[frame], zs[frame], marker = '.', color = 'y', s=1)





    
