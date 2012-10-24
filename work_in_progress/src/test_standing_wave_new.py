from wave_solver import solver
from numpy import *
import math as m
import os, glob, sys
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
    
C = 1
#Nx = 20
#Ny = 20
Lx = 4
Ly = 4
T = 4
c = 1.1
b = 0.05
w = pi
my = 2
mx = 2
try:
    version = sys.argv[1]
except:
    print "Please provide version (scalar or vectorized) on the command line."
    sys.exit(1)
q_const = 0.4
cx = mx*pi/Lx
cy = my*pi/Ly
dt = (1/sqrt(0.4))*(1/(sqrt(1/0.02**2 + 1/0.02**2)))


def exact_standing_wave(x,y,b,t,w,version="vectorized"):
    return m.exp(-b*t)*cos(mx*x*pi/Lx)*cos(my*y*pi/Ly)*m.cos(w*t)    
    

def test_standing_wave(version, w, Nx):
    Ny = Nx
   
    def q(x,y):
        return ones((len(x),len(y))) * q_const

    def f(x,y,t):
       return exp(-b*t)*cos(cx*x)*cos(cy*y)*(cos(w*t)*(q_const*cx**2 + q_const*cy**2 - w**2) + b*w*sin(w*t))

    def I(x,y):
        return cos(mx*x*pi/Lx)*cos(my*y*pi/Ly)

    def V(x,y):
        return -b*cos(mx*x*pi/Lx)*cos(my*y*pi/Ly)
        
    E_list, u, dx = solver(Lx, Ly, Nx, Ny, T, dt, c, I, q, V, f, b, version, w=w, exact=exact_standing_wave, make_plot=False)
    return E_list, dx
    
    
def compute_error(w, Nx):
    e_list, dx = test_standing_wave(version, w, Nx)
    E = sqrt(dx*dx*sum(e_list**2))
    return E, dx
       

Nx_list = [20, 40, 80, 160]
Error_list=[]
E_dx_list=[]
dx_list=[]
for nx in Nx_list:
    error, dx = compute_error(w, nx)
    dx_list.append(dx)
    Error_list.append(error)
    E_dx_list.append(error/dx**2)
    #print "Nx = %g, E = %g" % (nx, error)
m = len(Nx_list)
print "dx-list:"
print dx_list
print "E/dx**2: "
print E_dx_list
r = [log(Error_list[i-1]/Error_list[i])/log(dx_list[i-1]/dx_list[i]) for i in range(1, m, 1)]
print "Convergence rates: "
print r