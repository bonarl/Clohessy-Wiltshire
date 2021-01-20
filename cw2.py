#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Mar 23 13:38:54 2018

@author: bonar
@editor: aaronjohnsabu1999
"""

import math
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

plt.style.use('fivethirtyeight')

r0    = [0.0, 0.0, 0.0]       # initial starting position
rdot0 = [0.1, 0.0, 0.0]       # initial velocity (on release from ISS) (m/s)
R     =  6870 + 405
mu    =  398600.50
omega = math.sqrt(mu/R**3)

nframes01 = 1000
nframes02 =  200
nframes03 =  200
dt        =   20

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

def CW2(r0, rdot0, omega, t):
  x0    = r0[0]
  y0    = r0[1]
  z0    = r0[2]
  xdot0 = rdot0[0]
  ydot0 = rdot0[1]
  zdot0 = rdot0[2]
    
  A = [[         4.0 - 3.0*math.cos(omega*t), 0.0,                       0.0,           math.sin(omega*t)/omega,   2.0*(1.0-math.cos(omega*t))/omega,                     0.0],
       [ 6.0*math.sin(omega*t) - 6.0*omega*t, 1.0,                       0.0, 2.0*(math.cos(omega*t)-1.0)/omega, 4.0*math.sin(omega*t)/omega - 3.0*t,                     0.0],
       [                                 0.0, 0.0,         math.cos(omega*t),                               0.0,                                 0.0, math.sin(omega*t)/omega],
       [         3.0*math.sin(omega*t)*omega, 0.0,                       0.0,                 math.cos(omega*t),               2.0*math.sin(omega*t),                     0.0],
       [     6*omega*(math.cos(omega*t)-1.0), 0.0,                       0.0,           - 2.0*math.sin(omega*t),         4.0*math.cos(omega*t) - 3.0,                     0.0],
       [                                 0.0, 0.0, - omega*math.sin(omega*t),                               0.0,                                 0.0,       math.cos(omega*t)]]
  
  xt   = A[0][0]*r0[0] + A[0][1]*r0[1] + A[0][2]*r0[2] + A[0][3]*rdot0[0] + A[0][4]*rdot0[1] + A[0][5]*rdot0[2]
  yt   = A[1][0]*r0[0] + A[1][1]*r0[1] + A[1][2]*r0[2] + A[1][3]*rdot0[0] + A[1][4]*rdot0[1] + A[1][5]*rdot0[2]
  zt   = A[2][0]*r0[0] + A[2][1]*r0[1] + A[2][2]*r0[2] + A[2][3]*rdot0[0] + A[2][4]*rdot0[1] + A[2][5]*rdot0[2]
  xdot = A[3][0]*r0[0] + A[3][1]*r0[1] + A[3][2]*r0[2] + A[3][3]*rdot0[0] + A[3][4]*rdot0[1] + A[3][5]*rdot0[2]
  ydot = A[4][0]*r0[0] + A[4][1]*r0[1] + A[4][2]*r0[2] + A[4][3]*rdot0[0] + A[4][4]*rdot0[1] + A[4][5]*rdot0[2]
  zdot = A[5][0]*r0[0] + A[5][1]*r0[1] + A[5][2]*r0[2] + A[5][3]*rdot0[0] + A[5][4]*rdot0[1] + A[5][5]*rdot0[2]
  
  return([xt, yt, zt], [xdot, ydot, zdot])
        
xs = []
ys = []
zs = []
ds = []

d       =  0
d_max   = 50
stindex = nframes01

# Maneuver 1

for i in range(nframes01):
  t = dt*i
  if d < d_max:
    r_vec= CW2(r0, rdot0, omega, t)[0]
    rdot_vec = CW2(r0, rdot0, omega, t)[1]
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
    
# Maneuver 2
     
r0     = [x, y, z]
deltav = [-0.03, 0, 0.02]
for i in range(len(rdot)):
    rdot0[i] = rdot[i] + deltav[i]

for i in range(nframes02):
  t        = dt*i
  r_vec    = CW2(r0, rdot0, omega, t)[0]
  rdot_vec = CW2(r0, rdot0, omega, t)[1]
  x        = r_vec[0]
  y        = r_vec[1]
  z        = r_vec[2]
  xdot     = rdot_vec[0]
  ydot     = rdot_vec[1]
  zdot     = rdot_vec[2]
  rdot     = [xdot, ydot, zdot]
  d        = math.sqrt(x**2+y**2+z**2)
  ds.append(d)
  
  xs.append(x)
  ys.append(y)
  zs.append(z)
  
  v = math.sqrt(xdot**2+ydot**2+zdot**2)

# Maneuver 3

r0     = [x, y, z]
deltav = [0, -0.05, 0.03]
for i in range(len(rdot)):
  rdot0[i] = rdot[i] + deltav[i]  

for i in range(nframes03):
  t        = dt*i
  r_vec    = CW2(r0, rdot0, omega, t)[0]
  rdot_vec = CW2(r0, rdot0, omega, t)[1]
  x    = r_vec[0]
  y    = r_vec[1]
  z    = r_vec[2]
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
ax  = fig.add_subplot(111, projection = '3d', xlim = (xmin, xmax), ylim = (ymin, ymax), zlim = (zmin, zmax))
ax.set_xlabel('x')
ax.set_ylabel('y')
ax.set_zlabel('z')

frame = 0
while frame < stindex:
  frame +=1
  ax.scatter(xs[frame], ys[frame], zs[frame], marker = '.', color = 'g', alpha = 0.5, s=1)
  
while frame < nframes02+stindex:
  frame += 1
  ax.scatter(xs[frame], ys[frame], zs[frame], marker = '.', color = 'skyblue', s=1)
  
while frame < nframes03+nframes02+stindex:
  frame += 1
  ax.scatter(xs[frame], ys[frame], zs[frame], marker = '.', color = 'y', s=1)

plt.show()
