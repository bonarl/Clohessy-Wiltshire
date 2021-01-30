import math, plotly
import numpy as np
import matplotlib.pyplot as plt
import plotly.graph_objects as go
from mpl_toolkits.mplot3d import Axes3D

plt.style.use('fivethirtyeight')

r0   = [100.0, 50.0, 25.0]
rdot0 = [  1.0,  0.0,  0.0]
R     =  6870 + 405
mu    =  398600.50
omega = math.sqrt(mu/R**3)

nframes   = 2000
dt        =    2

def transformMatrix(omega, t):
  omegaT = omega*t
  cosOT  = math.cos(omegaT)
  sinOT  = math.sin(omegaT)
  
  A = [[        4.0 - 3.0*cosOT, 0.0,          0.0,           sinOT/omega,   2.0*(1.0-cosOT)/omega,         0.0],
       [ 6.0*sinOT - 6.0*omegaT, 1.0,          0.0, 2.0*(cosOT-1.0)/omega, 4.0*sinOT/omega - 3.0*t,         0.0],
       [                    0.0, 0.0,        cosOT,                   0.0,                     0.0, sinOT/omega],
       [        3.0*sinOT*omega, 0.0,          0.0,                 cosOT,               2.0*sinOT,         0.0],
       [    6*omega*(cosOT-1.0), 0.0,          0.0,           - 2.0*sinOT,         4.0*cosOT - 3.0,         0.0],
       [                    0.0, 0.0, -omega*sinOT,                   0.0,                     0.0,       cosOT]]
  return A

def coordAtT(r0, rdot0, omega, t):
  x0    = r0[0]
  y0    = r0[1]
  z0    = r0[2]
  xdot0 = rdot0[0]
  ydot0 = rdot0[1]
  zdot0 = rdot0[2]
    
  A = transformMatrix(omega, t)
  
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

# Rendezvous Maneuver
final_pos = [0.00, 0.00, 0.00]
T = dt*(nframes-1)
phi_r      = transformMatrix(omega, T)[0:3]
phi_r_r    = [arr[0:3] for arr in phi_r]
phi_r_rdot = [arr[3:6] for arr in phi_r]
deltav     = np.dot(np.linalg.inv(phi_r_rdot), (final_pos - np.dot(phi_r_r, r0))) - rdot0
print(deltav)
out = maneuver(nframes, r0, rdot0, deltav, omega)
r_vec      = out[0]
rdot_vec   = out[1]
xs.extend   (out[2])
ys.extend   (out[3])
zs.extend   (out[4])
ds.extend   (out[5])

xmin, xmax = min(xs), max(xs)
ymin, ymax = min(ys), max(ys)
zmin, zmax = min(zs), max(zs)
dmin, dmax = min(ds), max(ds)

marker_data = go.Scatter3d(x = xs[:nframes], y = ys[:nframes], z = zs[:nframes], mode = 'markers', marker = dict(size = 3))
fig = go.Figure(data = [marker_data])
plotly.offline.plot(fig)
