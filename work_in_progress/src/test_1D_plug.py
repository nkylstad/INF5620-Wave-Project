from wave_solver import solver
from numpy import *

C = 1
Nx = 50
Ny = 50
Lx = 4
Ly = 4
T = 10
c = 1
b = 0.0
dt = -1.0
sigma = Lx/10.
#xc = Lx/2.
xc = 0

def test_1D_plug(version):
    
    def q(x,y):
        return ones((len(x),len(y)))

    def V(x,y):
        return zeros((len(x), len(y)))
    
    if version == "scalar":
        def I(x,y):
            return 0 if abs(x-xc) > sigma else 2
    else: 
        def I(xv, yv):
            Iv = zeros((len(xv), len(yv)))
            for i in range(len(xv)):
                for j in range(len(yv)):
                    if abs(xv[j, i] - xc) > sigma:
                        Iv[i, j] = 0
                    else:
                        Iv[i, j] = 2
                   #print "iv:",Iv[i, j]
                        
            return Iv
        
    def f(x,y,t):
        return zeros((len(x), len(y)))
        
    E, u, dx = solver(Lx, Ly, Nx, Ny, T, dt, c, I, q, V, f, b, version, oneD=True)
    

#test_constant_1D("scalar")    
test_1D_plug("vectorized")