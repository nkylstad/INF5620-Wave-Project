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
dt = (1/sqrt(Lx*Ly))*(1/(sqrt(1/0.02**2 + 1/0.02**2)))
w = pi
my = 2
mx = 2
h = 0.1
try:
    version = sys.argv[1]
except:
    print "Please provide version (scalar or vectorized) on the command line."
    sys.exit(1)
    
cx = mx*pi/Lx
cy = my*pi/Ly


def exact_manufactured_solution(x,y,b,t,w,version="vectorized"):
    return exp(-b*t)*cos(mx*x*pi/Lx)*cos(my*y*pi/Ly)*cos(w*t)
     
          
def test_manufactured_solution(version, w, Nx):
   
    Ny = Nx

    def q(x,y):
        return x*y
    
    def f(x,y,t):
        f_val = exp(-b*t)*(cos(cx*x)*cos(cy*y)*(-w**2*cos(w*t) + b*w*sin(w*t)) + \
            cx*y*cos(cy*y)*cos(w*t)*(cx*x*cos(cx*x) + sin(cy*y)) + \
            cy*x*cos(cx*x)*cos(w*t)*(cy*y*cos(cy*y) + sin(cx*x)))
        return f_val
        
    def I(x,y):
        if version=="scalar":   
            return m.cos(cx*x)*m.cos(cy*y)
        else:
            return cos(cx*x)*cos(cy*y)
        
    def V(x,y):
        return -b*cos(cx*x)*cos(cy*y)
        
    E_list, u, dx = solver(Lx, Ly, Nx, Ny, T, dt, c, I, q, V, f, b, version, w=w, exact=exact_manufactured_solution, make_plot=False)
    return E_list, dx
    
    
def compute_error(w, Nx):
    e_list, dx = test_manufactured_solution(version, w, Nx)
    E = sqrt(dx*dx*sum(e_list**2))
    return E, dx
       

Nx_list = [20, 40, 80]
Error_list=[]
dx_list = []
E_dx_list = []
for nx in Nx_list:
    error, dx = compute_error(w, nx)
    Error_list.append(error)
    dx_list.append(dx)
    E_dx_list.append(error/(dx**2))
m = len(Nx_list)
print "dx-list:"
print dx_list
print "E/dx**2: "
print E_dx_list
print "Error_list:"
print Error_list
r = [log(Error_list[i-1]/Error_list[i])/log(dx_list[i-1]/dx_list[i]) for i in range(1, m, 1)]
print "Convergence rates: "
print r