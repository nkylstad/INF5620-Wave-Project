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
    
    def q(x,y):
        return ones((len(x),len(y)))*0.8
    
    def V(x,y):
        if version == "scalar":
            return 0
        else:
            return zeros((len(x),len(y)))
        
        
    def I(x,y):
        return 1.2
        
    
    def f(x,y,t):
        if version == "scalar":
            return 0
        else:
            return zeros((len(x),len(y)))
        
    
    E, u, dx = solver(Lx,Ly,Nx,Ny,T,dt,c,I,q,V,f,b,version)
    if version == "scalar":
        print u
    else:
        print u[1:-1,1:-1]
 
print "Scalar:"
test_constant("scalar")
print "Vectorized:"
test_constant("vectorized")