from wave_solver import solver
from numpy import *

C = 1
Nx = 50
Ny = 50
Lx = 10
Ly = 10
T = 10
c = 1
b = 0.1
dt = -1.0
sigma = Lx/10.
xc = Lx/2.
#xc = 0

def test_1D_plug(version):
    
    q = ones((Nx+1,Ny+1))
    
    def V(x,y):
        if version=="scalar":
            return 0
        else: 
            return zeros((Nx+1, Ny+1))
    
    if version == "scalar":
        def I(x,y):
            return 0 if abs(x-xc) > sigma else 2
    else: 
        def I(xv, yv):
            Iv = zeros((Nx+1, Ny+1))
            for i in range(Nx+1):
                for j in range(Ny+1):
                    if abs(xv[j, i] - xc) > sigma:
                        Iv[i, j] = 0
                    else:
                        Iv[i, j] = 2
                    #print "iv:",Iv[i, j]
                        
            return Iv
        
    def f(x,y,t):
        if version=="scalar":
            return 0
        else:
            return zeros((Nx+1, Ny+1))
        
    u = solver(Lx, Ly, Nx, Ny, T, dt, c, I, q, V, f, b, version, oneD=True)
    

#test_constant_1D("scalar")    
test_1D_plug("vectorized")