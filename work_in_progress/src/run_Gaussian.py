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
T = 1
Lx = 4
Ly = 4
c = 1.1
b = 0.1
sigma = 0.8
factor = 0
Ba = 20
Bsx = 0.25
Bsy = 0.25
g = 9.81
H0=50

def run_Gaussian(version):
    def I(x,y):
        """Gaussian peak at (0, Ly/2)."""
        # return 2*exp(-0.5*((x-Lx/2.0)/(sigma))**2 - 0.5*((y-Ly/2.0)/(sigma))**2)
        return 2*exp(-0.5*((x)/(sigma))**2 - 0.5*((y-Ly/2.0)/(sigma))**2)

    def B(x,y):
        return Ba*exp(-((x-Lx/2)/Bsx)**2 - ((y-Ly/2)/Bsy)**2)

    def V(x,y):
        return ones((len(x),len(y)))*factor
            
    
    def q(x,y):
        return g*(H0-B(x,y))

    # def q(x, y):
    #     q = ones((len(x),len(y)))
    #     return 0.3*q
    
    def f(x,y,t):
        return zeros((len(x),len(y)))
    
    E, u, dx = solver(Lx, Ly, Nx, Ny, T, dt, c, I, q, V, f, b, version, make_plot=True)
    
run_Gaussian("vectorized")
