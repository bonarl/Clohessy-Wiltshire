import math
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

plt.style.use('fivethirtyeight')

r0    = [0.0, 0.0, 0.0]       # initial starting position
rdot0 = [0.1, 0.0, 0.0]       # initial velocity (on release from ISS) (m/s)
R     =  6870 + 405
mu    =  398600.50
omega = math.sqrt(mu/R**3)

nframes01 =  200
nframes02 =  200
nframes03 =  200
dt        =   20

def coordAtT(r0, rdot0, omega, t):
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
  
  return [xt, yt, zt], [xdot, ydot, zdot]

def maneuver(frames, r0, rdot0, deltav, omega):
  xs, ys, zs, ds = [], [], [], []
  for i, vectorElem in enumerate(rdot0):
    rdot0[i] = vectorElem + deltav[i]
  for i in range(frames):
    t        = dt*i
    r, rdot  = coordAtT(r0, rdot0, omega, t)
    xs.append(r[0])
    ys.append(r[1])
    zs.append(r[2])
    ds.append(np.linalg.norm(r))
    v = math.sqrt(np.linalg.norm(rdot))
  return r, rdot, xs, ys, zs, ds, v
  
xs = []
ys = []
zs = []
ds = []

# Maneuver 1
out = maneuver(nframes01, r0,    rdot0,    [ 0.00,  0.00, 0.00], omega)
r_vec    = out[0]
rdot_vec = out[1]
xs.extend (out[2])
ys.extend (out[3])
zs.extend (out[4])
ds.extend (out[5])

# Maneuver 2
out = maneuver(nframes02, r_vec, rdot_vec, [-0.03,  0.00, 0.02], omega)
r_vec    = out[0]
rdot_vec = out[1]
xs.extend (out[2])
ys.extend (out[3])
zs.extend (out[4])
ds.extend (out[5])

# Maneuver 3
out = maneuver(nframes03, r_vec, rdot_vec, [ 0.00, -0.05, 0.03], omega)
r_vec    = out[0]
rdot_vec = out[1]
xs.extend (out[2])
ys.extend (out[3])
zs.extend (out[4])
ds.extend (out[5])

xmin, xmax = min(xs), max(xs)
ymin, ymax = min(ys), max(ys)
zmin, zmax = min(zs), max(zs)
dmin, dmax = min(ds), max(ds)

fig = plt.figure()
ax  = fig.add_subplot(111, projection = '3d', xlim = (xmin, xmax), ylim = (ymin, ymax), zlim = (zmin, zmax))
ax.set_xlabel('x')
ax.set_ylabel('y')
ax.set_zlabel('z')

for frame in range(nframes01 + nframes02 + nframes03):
  if frame < nframes01:
    ax.scatter(xs[frame], ys[frame], zs[frame], marker = '.', color = 'g', s=1)
  elif frame < nframes01 + nframes02:
    ax.scatter(xs[frame], ys[frame], zs[frame], marker = '.', color = 'b', s=1)
  else:
    ax.scatter(xs[frame], ys[frame], zs[frame], marker = '.', color = 'r', s=1)

plt.show()