from numpy import *
import math as m
from mayavi import mlab
from scitools.easyviz import *

c = 1.1
T = 2
Nx = 20
Ny = 20
Lx = 4
Ly = 4
b = 0
w = pi
mx = 2
my = 2
sigma = 0.8
x = linspace(0, Lx, Nx+1)
y = linspace(0, Ly, Ny+1)
X,Y = meshgrid(x,y)
dx = x[1] - x[0]
dy = y[1] - y[0]
dt = (1/c)*(1/sqrt(1/dx**2 + 1/dy**2))
N = int(round(float(T/dt)))
t = linspace(0, T, N+1)

def I(x,y):
    """Gaussian peak at (Lx/2, Ly/2)."""
    return 2*exp(-0.5*((x-Lx/2.0)/(sigma))**2 - 0.5*((y-Ly/2.0)/(sigma))**2)
    
q = I(X,Y) - 5


def exact_standing_wave(x,y,b,t,w,version="vectorized"):
    if version=="scalar":
        return m.exp(-b*t)*m.cos(mx*x*pi/Lx)*m.cos(my*y*pi/Ly)*m.cos(w*t)
    else:
        return m.exp(-b*t)*cos(mx*x*pi/Lx)*cos(my*y*pi/Ly)*m.cos(w*t)
        
q_axes = (0, Lx+2, 0, Ly+2, -5, -1.5)
u0 = exact_standing_wave(X, Y, b, t[0], w)
s = mlab.mesh(X, Y, u0,colormap='cool', opacity=0.5)  # color=(0.1,0.7,0.8)
qp = mlab.mesh(X,Y,q, color=(0.2,0.2,0.2), opacity=0.8, reset_zoom=False)
mlab.options.offscreen = True
mlab.savefig("u%.3d.png" % 0)
mlab.clf()

for n in range(1,N+1):
    u = exact_standing_wave(X, Y, b, t[n], w)
    s = mlab.mesh(X,Y,u,colormap='cool', opacity=0.8)
    qp = mlab.mesh(X,Y,q, color=(0.7,0.7,0.7), opacity=0.7, reset_zoom=False)
    mlab.options.offscreen = True
    mlab.savefig("u%.3d.png" % n)
    mlab.clf()
    
import os, glob

movie("u*.png")
for i in glob.glob("u*.png"):
    os.remove(i)
