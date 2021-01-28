import math, plotly
import numpy as np
import matplotlib.pyplot as plt
import plotly.graph_objects as go
from mpl_toolkits.mplot3d import Axes3D

plt.style.use('fivethirtyeight')

r_0   = [100.0, 50.0, 25.0]
rdot0 = [  1.0,  0.0,  0.0]
R     =  6870 + 405
mu    =  398600.50
omega = math.sqrt(mu/R**3)

nframes =  2000
dt      =  4000.0/nframes

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
  
  xt   = A[0][0]*x0 + A[0][1]*y0 + A[0][2]*z0 + A[0][3]*xdot0 + A[0][4]*ydot0 + A[0][5]*zdot0
  yt   = A[1][0]*x0 + A[1][1]*y0 + A[1][2]*z0 + A[1][3]*xdot0 + A[1][4]*ydot0 + A[1][5]*zdot0
  zt   = A[2][0]*x0 + A[2][1]*y0 + A[2][2]*z0 + A[2][3]*xdot0 + A[2][4]*ydot0 + A[2][5]*zdot0
  xdot = A[3][0]*x0 + A[3][1]*y0 + A[3][2]*z0 + A[3][3]*xdot0 + A[3][4]*ydot0 + A[3][5]*zdot0
  ydot = A[4][0]*x0 + A[4][1]*y0 + A[4][2]*z0 + A[4][3]*xdot0 + A[4][4]*ydot0 + A[4][5]*zdot0
  zdot = A[5][0]*x0 + A[5][1]*y0 + A[5][2]*z0 + A[5][3]*xdot0 + A[5][4]*ydot0 + A[5][5]*zdot0
  
  return [xt, yt, zt], [xdot, ydot, zdot]

def maneuver(frames, r0, rdot0, deltav, omega):
  xs, ys, zs, ds = [], [], [], []
  rdot0 = np.add(rdot0, deltav)
  for i in range(frames):
    t        = dt*i
    r, rdot  = coordAtT(r0, rdot0, omega, t)
    xs.append(r[0])
    ys.append(r[1])
    zs.append(r[2])
    ds.append(np.linalg.norm(r))
    v = math.sqrt(np.linalg.norm(rdot))
  return r, rdot, xs, ys, zs, ds, v

# Glideslope Rendezvous Maneuver
def glideslope(N):
  xs = []
  ys = []
  zs = []
  ds = []

  r_T       = [0.00, 0.00, 0.00]
  T         = dt*(nframes-1)
  deltaT    = T/N
  
  u_rho     = [(r_0[i] - r_T[i]) for i in range(len(r_T))]
  rho_0     = np.linalg.norm(u_rho)
  u_rho    /= rho_0
  r_vec     = r_0
  rdot_vec  = rdot0
  
  deltavs   = []
  energy    = 0.0
  
  for firing in range(N):
    r_m   = r_vec
    rho_m = np.linalg.norm([(r_m[i] - r_T[i]) for i in range(len(r_T))])
    r_m1  = [(r_T[i] + rho_0*(N-firing-1)*u_rho[i]/N) for i in range(len(r_T))]
    
    phi_r      = transformMatrix(omega, deltaT)[0:3]
    phi_r_r    = [arr[0:3] for arr in phi_r]
    phi_r_rdot = [arr[3:6] for arr in phi_r]
    deltav     = np.dot(np.linalg.inv(phi_r_rdot), (r_m1 - np.dot(phi_r_r, r_m))) - rdot_vec
    deltavs.append(np.linalg.norm(deltav))
    energy    += deltavs[-1]
    #energy    += np.linalg.norm(deltav)
    
    out        = maneuver(int(nframes/N), r_vec, rdot_vec, deltav, omega)
    r_vec      = out[0]
    rdot_vec   = out[1]
    xs.extend   (out[2])
    ys.extend   (out[3])
    zs.extend   (out[4])
    ds.extend   (out[5])
  
  return xs, ys, zs, ds, deltavs, energy

def simulate(xs, ys, zs, ds):
  xmin, xmax = min(xs), max(xs)
  ymin, ymax = min(ys), max(ys)
  zmin, zmax = min(zs), max(zs)
  dmin, dmax = min(ds), max(ds)
  
  marker_data = go.Scatter3d(x = xs[:nframes], y = ys[:nframes], z = zs[:nframes], mode = 'markers', marker = dict(size = 3))
  fig = go.Figure(data = [marker_data])
  plotly.offline.plot(fig)

energies = []
minN  =  1
maxN  = 10
allNs = []
for N in range(minN,maxN+1):
  #if nframes%N==0 and nframes > N:
  if True:
    xs, ys, zs, ds, deltavs, energy = glideslope(N)
    print(energy, deltavs)
    print(N)
    allNs.append(N)
    energies.append(energy)
    simulate(xs, ys, zs, ds)

plt.scatter(allNs, energies)
plt.xscale('log')
plt.show()