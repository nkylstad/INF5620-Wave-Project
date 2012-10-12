from numpy import *
from wave_solver import solver
import math as m
import os, glob

#os.remove("initial.txt")
#os.remove("u0.txt")
#for i in glob.glob("texttmp*.txt"):
    #os.remove(i)
#for j in glob.glob("wtmp*.png"):
    #os.remove(j)


dt = -1
Nx = 20
Ny = 20
T = 10
Lx = 4
Ly = 4
c = 1.1
b = 0.1
sigma = 0.8

def run_Gaussian(version):
    def I(x,y):
        """Gaussian peak at (Lx/2, Ly/2)."""
        return 2*exp(-0.5*((x-Lx/2.0)/(sigma))**2 - 0.5*((y-Ly/2.0)/(sigma))**2)
    version = version
    V = ones((Nx+1,Ny+1))
    V = V*0.01
    q = ones((Nx+1,Ny+1))
    q = 0.3*q
    
    if version == "scalar":
        def f(x,y,t):
            return 0
    else:
        def f(x,y,t):
            return zeros((Nx+1,Ny+1))
    
    
    E = solver(Lx, Ly, Nx, Ny, T, dt, c, I, q, V, f, b, version, make_plot=True)
    
run_Gaussian("vectorized")
