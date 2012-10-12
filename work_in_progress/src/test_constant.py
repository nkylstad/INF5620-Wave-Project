from wave_solver import solver
from numpy import zeros, ones

dt = -1
Nx = 50
Ny = 50
T = 10
Lx = 10
Ly = 10
c = 1.1
b = 0.1
sigma = 0.8

def test_constant(version):
    q = ones((Nx+1,Ny+1))*0.8
    
    if version == "scalar":
        def V(x,y):
            return 0
    else:
        V = zeros((Nx+1,Ny+1))
        
        
    def I(x,y):
        return 1.2
        
    if version == "scalar":
        def f(x,y,t):
            return 0
    else:
        f = zeros((Nx+1,Ny+1))
        
    
    E = solver(Lx,Ly,Nx,Ny,T,dt,c,I,q,V,f,b,version)
 
#test_constant("scalar")
#test_constant("vectorized")