from numpy import *
import sys
from scitools.easyviz import *
import time
#from matplotlib.pyplot import *


def solver(I, V, q, f, Lx, Ly, Nx, Ny, c, dt, T, b, version, make_plot=True):
    
    x = linspace(0,Lx,Nx+1)  # mesh points in x-direction
    y = linspace(0,Ly,Ny+1)  # mesh points in y-direction
    xv = x.reshape((x.size, 1))
    yv = y.reshape((x.size, 1))
    dx = x[0] - x[1]
    dy = y[0] - y[1]
    
    stability_limit = (1/float(c))*(1/sqrt(1/dx**2 + 1/dy**2))
    if dt <= 0:
        dt = 1*stability_limit
    elif dt > stability_limit:
        print 'error: dt=%g exceeds the stability limit %g' % \
              (dt, stability_limit)


    N = int(round(T/float(dt)))
    t = linspace(0,N*dt, N+1)  # mesh points in time
    Cx2 = (dt/dx)**2; Cy2 = (dt/dy)**2
    dt2 = dt**2
    c1 = 1/(1 + ((b*dt)/2))
    c2 = 1 - (b*dt)/2
    
    up = zeros((Nx+1,Ny+1))
    u = zeros((Nx+1, Ny+1))
    um = zeros((Nx+1, Ny+1))
    t0 = time.clock()
    # Initial conditions
    if version == "scalar":
        for i in range(Nx+1):
            for j in range(Ny+1):
                u[i,j] = I(x[i],y[j])
    else:  # use vectorized version
        u[:,:] = I(xv, yv)
            
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
            
            up[i][j] = c1*(2*u[i][j] - c2*(u[i][j] - dt*V(x[i],y[j])) + 0.5*Cx2*((q(x[i],y[j]) + q(x[spx],y[j]))*\
            (u[spx][j] - u[i][j]) - (q(x[smx],y[j]) + q(x[i],y[j]))*(u[i][j] - u[smx][j])) + \
            0.5*Cy2*((q(x[i],y[j]) + q(x[i],y[spy]))*(u[i][spy] - u[i][j]) - (q(x[i],y[j]) + q(x[i],y[smy]))*\
            (u[i][j] - u[i][smy])) + dt2*f(x[i],y[j],t[0])[i,j])
            
    # update arrays one timestep:
    um = u
    u = up
    for n in range(1, N+1):
        if version == "scalar":
            up, u, um = advance_scalar(Nx, Ny, x, y, t[n], up, u, um, c1, c2, Cx2, Cy2, dt2, q, f)
        else:
            qv = q(xv,yv)
            fv = f(xv,yv,t)
            up, u, um = advance_vectorized(Nx, Ny, x, y, up, u, um, c1, c2, Cx2, Cy2, dt2, qv, fv)
        if make_plot:
            plot_3D(x,y,u)
    t1 = time.clock()
    print "Finished! Used %g seconds" % (t1-t0)
        
    return True
    
def advance_scalar(Nx, Ny, x, y, t, up, u, um, c1, c2, Cx2, Cy2, dt2, q, f):
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
            up[i][j] = c1*(2*u[i][j] - c2*um[i][j] + 0.5*Cx2*((q(x[i],y[j]) + q(x[spx],y[j]))*\
            (u[spx][j] - u[i][j]) - (q(x[smx],y[j]) + q(x[i],y[j]))*(u[i][j] - u[smx][j])) + \
            0.5*Cy2*((q(x[i],y[j]) + q(x[i],y[spy]))*(u[i][spy] - u[i][j]) - (q(x[i],y[j]) + q(x[i],y[smy]))*\
            (u[i][j] - u[i][smy])) + dt2*f(x[i],y[j],t))

    um = u
    u = up
    return up, u, um

def advance_vectorized(Nx, Ny, x, y, up, u, um, c1, c2, Cx2, Cy2, dt2, q, f):
    
    up[1:-1,1:-1] = c1*(2*u[1:-1,1:-1] - c2*um[1:-1,1:-1] + \
        0.5*Cx2*((q[:-2] + q[1:-1])*(u[:-2,1:-1] - 2*um[1:-1,1:-1]) -\
         (q[1:-1] + q[2:])*(u[1:-1,1:-1] - 2*um[2:,1:-1])) + \
         0.5*Cy2*((q[:-2] + q[1:-1])*(u[1:-1,:-2] - 2*um[1:-1,1:-1]) -\
         (q[1:-1] + q[2:])*(u[1:-1,1:-1] - 2*um[1:-1,2:])) + \
         dt2*f[1:-1,1:-1])
    
    #Boundary conditions:
    up[0,:] = 0
    up[Nx,:] = 0
    up[:,Ny] = 0
    up[:,0] = 0
    
    um = u
    u = up
    return um, u, up

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
    version = "vectorized"
    #version = "scalar"
    def I(x,y):
        """Gaussian peak at (Lx/2, Ly/2)."""
        return exp(-0.5*(x-Lx/2.0)**2 - 0.5*(y-Ly/2.0)**2)
    
    def V(x,y):
        return 0
    
    def q(x,y):
        return 1*x + 1*y
    
    if version == "scalar":
        def f(x,y,t):
            return 0
    else:
        def f(x,y,t):
            return zeros((Nx+1,Ny+1))
    
    solver(I, V, q, f, Lx, Ly, Nx, Ny, c, dt, T, b, version, make_plot=True)
    

run_Gaussian()
print "No syntax errors. Yeey!"

   