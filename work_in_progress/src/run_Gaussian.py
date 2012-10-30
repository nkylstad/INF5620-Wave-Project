from numpy import *
from wave_solver import solver
import math as m
import os, glob

dt = -1
Nx = 100
Ny = 100
T = 3
Lx = 40
Ly = 40
c = 1.1
b = 0.1
sigma = 4.0
factor = 0
height = 6
Ba = 20
Bsx = 6.0
Bsy = 6.0
g = 9.81
H0=25

f = open("H0.txt", 'w')
f.write("%g" % H0)
f.close()


def run_Gaussian(version):
    def I(x,y):
        """Gaussian peak at (0, Ly/2)."""
        # return 2*exp(-0.5*((x-Lx/2.0)/(sigma))**2 - 0.5*((y-Ly/2.0)/(sigma))**2)
        return height*exp(-0.5*((x)/(sigma))**2 - 0.5*((y-Ly/2.0)/(sigma))**2)

    def B(x,y):
        return Ba*exp(-((x-Lx/2)/Bsx)**2 - ((y-Ly/2)/Bsy)**2)
        #return -x

    def V(x,y):
        return ones((len(x),len(y)))*factor
            
    
    def q(x,y):
        return g*(H0-B(x,y))

    def f(x,y,t):
        return zeros((len(x),len(y)))
    
    E, u, dx = solver(Lx, Ly, Nx, Ny, T, dt, c, I, q, V, f, b, version, B=B, make_plot=True)
    
#run_Gaussian("vectorized")

def run_Gaussian_flathill(version):
    def I(x,y):
        """Gaussian peak at (0, Ly/2)."""
        # return 2*exp(-0.5*((x-Lx/2.0)/(sigma))**2 - 0.5*((y-Ly/2.0)/(sigma))**2)
        return height*exp(-0.5*((x)/(sigma))**2 - 0.5*((y-Ly/2.0)/(sigma))**2)

    def B(x,y):
        return x-H0

    def V(x,y):
        return ones((len(x),len(y)))*factor
            
    
    def q(x,y):
        return g*(H0-B(x,y))

    def f(x,y,t):
        return zeros((len(x),len(y)))
    
    E, u, dx = solver(Lx, Ly, Nx, Ny, T, dt, c, I, q, V, f, b, version, B=B, make_plot=True)

run_Gaussian_flathill("vectorized")