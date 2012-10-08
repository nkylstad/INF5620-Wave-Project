from numpy import linspace, zeros
from math import *
import sys
from scitools.easyviz import *
#from matplotlib.pyplot import *


def solver(I, V, q, f, Lx, Ly, Nx, Ny, c, dt, T, b, make_plot=True):
    
    x = linspace(0,Lx,Nx+1)  # mesh points in x-direction
    y = linspace(0,Ly,Ny+1)  # mesh points in y-direction
    dx = x[0] - x[1]
    dy = y[0] - y[1]
    
    stability_limit = (1/float(c))*(1/sqrt(1/dx**2 + 1/dy**2))
    if dt <= 0:
        dt = 1*stability_limit
    elif dt > stability_limit:
        print 'error: dt=%g exceeds the stability limit %g' % \
              (dt, stability_limit)


    N = int(round((T/float(dt))))
    t = linspace(0,N*dt, N+1)  # mesh points in time
    Cx2 = (dt/(sqrt(2)*dx))**2; Cy2 = (dt/(sqrt(2)*dy))**2
    c1 = 1/(1 + ((b*dt)/2))
    c2 = 1 + (b*dt)/2
    
    up = zeros((Nx+1,Ny+1))
    u = zeros((Nx+1, Ny+1))
    um = zeros((Nx+1, Ny+1))
    
    # Initial conditions
    for i in range(Nx+1):
        for j in range(Ny+1):
            u[i][j] = I(x[i],y[j])
            
    # First step:
    for i in range(Nx+1):
        for j in range(Ny+1):
            #Boundary conditions:
            if i == 0: smx = i+1; 
            else: smx = i-1
            if i == Nx: spx = i-1; 
            else: spx = i+1
            if j == 0: smy = j+1; 
            else: smy = j-1
            if j == Ny: spy = j-1; 
            else: spy = j+1
            
            up[i][j] = c1*(2*u[i][j] - c2*(u[i][j] - dt*V(x[i],y[j])) + Cx2*((q(x[i],y[j]) + q(x[spx],y[j]))*\
            (u[spx][j] - u[i][j]) - (q(x[smx],y[j]) + q(x[i],y[j]))*(u[i][j] - u[smx][j])) + \
            Cy2*((q(x[i],y[j]) + q(x[i],y[spy]))*(u[i][spy] - u[i][j]) - (q(x[i],y[j]) + q(x[i],y[smy]))*\
            (u[i][j] - u[i][smy])) + f(x[i],y[j],t[0]))
            
    # update arrays one timestep:
    um = u
    u = up
    for n in range(1, N+1):
        up, u, um = advance_scalar(Nx, Ny, x, y, t[n], up, u, um, c1, c2, Cx2, Cy2, q, f)
        if make_plot:
            plot_3D(x, y, u)
        
    return True
    
def advance_scalar(Nx, Ny, x, y, t, up, u, um, c1, c2, Cx2, Cy2, q, f):
    for i in range(Nx+1):
        for j in range(Ny+1):
            #Boundary conditions:
            if i == 0: smx = i+1; 
            else: smx = i-1
            if i == Nx: spx = i-1; 
            else: spx = i+1
            if j == 0: smy = j+1; 
            else: smy = j-1
            if j == Ny: spy = j-1; 
            else: spy = j+1
            
            up[i][j] = c1*(2*u[i][j] - c2*um[i][j] + Cx2*((q(x[i],y[j]) + q(x[spx],y[j]))*\
            (u[spx][j] - u[i][j]) - (q(x[smx],y[j]) + q(x[i],y[j]))*(u[i][j] - u[smx][j])) + \
            Cy2*((q(x[i],y[j]) + q(x[i],y[spy]))*(u[i][spy] - u[i][j]) - (q(x[i],y[j]) + q(x[i],y[smy]))*\
            (u[i][j] - u[i][smy])) + f(x[i],y[j],t))
            
    um = u
    u = up
            
    return up, u, um

def plot_3D(x, y, u):
    surfc(x,y,u)

def run_Gaussian():
    dt = -1
    Nx = 100
    Ny = 100
    T = 10
    Lx = 10
    Ly = 10
    c = 1.0
    b = 1.0
    def I(x,y):
        """Gaussian peak at (Lx/2, Ly/2)."""
        return exp(-0.5*(x-Lx/2.0)**2 - 0.5*(y-Ly/2.0)**2)
    
    def V(x,y):
        return 0
    
    def q(x,y):
        return 1
    
    def f(x,y,t):
        return 0
    
    solver(I, V, q, f, Lx, Ly, Nx, Ny, c, dt, T, b, make_plot=True)
    

run_Gaussian()
print "No syntax errors. Yeey!"

   