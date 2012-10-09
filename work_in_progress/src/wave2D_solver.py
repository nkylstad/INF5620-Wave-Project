from numpy import *
import sys
#from scitools.easyviz import movie
import time
from mayavi import mlab
#from matplotlib.pyplot import *


def solver(I, V, q, f, Lx, Ly, Nx, Ny, c, dt, T, b, version, make_plot=True):
    
    x = linspace(0,Lx,Nx+1)  # mesh points in x-direction
    y = linspace(0,Ly,Ny+1)  # mesh points in y-direction
    xv = x[:,newaxis] # for vectorized function evaluations
    yv = y[newaxis,:]
    X,Y = meshgrid(x,y)
    dx = x[1] - x[0]
    dy = y[1] - y[0]
    
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
                um[i,j] = I(x[i],y[j])
    else:  # use vectorized version
        um[:,:] = I(X, Y)
            
    # First step:
    #up = u.copy()
    
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
            
            
            u[i,j] = um[i,j] + dt*c2*V[i,j] + 0.5*Cx2*((q[i,j] + q[spx,j])*\
            (um[spx,j] - um[i,j]) - (q[smx,j] + q[i,j])*(um[i,j] - um[smx,j])) + \
            0.5*Cy2*((q[i,j] + q[i,spy])*(um[i,spy] - um[i,j]) - (q[i,j] + q[i,smy])*\
            (um[i,j] - um[i,smy])) #+ 0.5*dt2*f(x[i],y[j],t[0])[i,j])
            
            
    # update arrays one timestep:
    #u = um.copy()
    #u = up.copy()
    s = mlab.mesh(X[1:-1,1:-1], Y[1:-1,1:-1], u[1:-1,1:-1])
    
    for n in range(2, N+1):
        if version == "scalar":
            u, um = advance_scalar(Nx, Ny, x, y, t[n], u, um, c1, c2, Cx2, Cy2, dt2, q, f)
            if make_plot:
		plot_3D(x,y,u,s)
        else:
            u, um = advance_vectorized(Nx, Ny, x, y, u, um, c1, c2, Cx2, Cy2, dt2, q, f)
	    if make_plot:
		plot_3D(X,Y,u,s)
    t1 = time.clock()
    print "Finished! Used %g seconds" % (t1-t0)
        
    return dt
    
def advance_scalar(Nx, Ny, x, y, t, u, um, c1, c2, Cx2, Cy2, dt2, q, f):
    up = u.copy()
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
            up[i,j] = c1*(2*u[i,j] - c2*um[i,j] + 0.5*Cx2*((q[i,j] + q[spx,j])*\
            (u[spx,j] - u[i,j]) - (q[smx,j] + q[i,j])*(u[i,j] - u[smx,j])) + \
            0.5*Cy2*((q[i,j] + q[i,spy])*(u[i,spy] - u[i,j]) - (q[i,j] + q[i,smy])*\
            (u[i,j] - u[i,smy]))) #+ dt2*f(x[i],y[j],t))

    um = u.copy()
    u = up.copy()
    return u, um

def advance_vectorized(Nx, Ny, x, y, u, um, c1, c2, Cx2, Cy2, dt2, q, f):
    up = u.copy()
    up[1:-1,1:-1] = c1*(2*u[1:-1,1:-1] - c2*um[1:-1,1:-1] + \
	  0.5*Cx2*((q[1:-1,1:-1] + q[2:,1:-1])*(u[2:,1:-1]-u[1:-1,1:-1]) - \
	  (q[1:-1,1:-1] + q[:-2,1:-1])*(u[1:-1,1:-1] - u[:-2,1:-1])) + \
	  0.5*Cy2*((q[1:-1,1:-1] + q[1:-1,2:])*(u[1:-1,2:] - u[1:-1,1:-1]) - \
	  (q[1:-1,1:-1] + q[1:-1,:-2])*(u[1:-1,1:-1] - u[1:-1,:-2])))
    
    """
    up[1:-1,1:-1] = 2*u[1:-1,1:-1] - um[1:-1,1:-1] + \
          Cx2*((u[2:,1:-1] - 2*u[1:-1,1:-1] + u[:-2,1:-1])) + \
          Cy2*((u[1:-1,2:] - 2*u[1:-1,1:-1] + u[1:-1,:-2]))
    """
    
    #boundary conditions:
    """
    u[0,:] = u[1,:]
    u[-1,:] = u[-2,:]
    u[:,0] = u[:,1]
    u[:,-1] = u[:,-2]
  
    """
    # x = 0
    up[0,1:-1] = c1*(2*u[0,1:-1] - c2*um[0,1:-1] + \
	  0.5*Cx2*((q[0,1:-1] + q[1,1:-1])*(u[1,1:-1] - u[0,1:-1]) - \
	  (q[0,1:-1] + q[1,1:-1])*(u[0,1:-1] - u[1,1:-1])) + \
	  0.5*Cy2*((q[0,1:-1] + q[0,2:])*(u[0,2:] - u[0,1:-1]) - \
	  (q[0,1:-1] + q[0,:-2])*(u[0,1:-1] - u[0,:-2])))
    
    # x = Lx
    up[-1,1:-1] = c1*(2*u[-1,1:-1] - c2*um[-1,1:-1] + \
	  0.5*Cx2*((q[-1,1:-1] + q[-2,1:-1])*(u[-2,1:-1] - u[-1,1:-1]) - \
	  (q[-1,1:-1] + q[-2,1:-1])*(u[-1,1:-1] - u[-2,1:-1])) + \
	  0.5*Cy2*((q[-1,1:-1] + q[-1,2:])*(u[-1,2:] - u[-1,1:-1]) - \
	  (q[-1,1:-1] + q[-1,:-2])*(u[-1,1:-1] - u[-1,:-2])))

    #y = 0
    up[1:-1,0] = c1*(2*u[1:-1,0] - c2*um[1:-1,0] + \
	  0.5*Cx2*((q[1:-1,0] + q[2:,0])*(u[2:,0] - u[1:-1,0]) - \
	  (q[1:-1,0] + q[:-2,0])*(u[1:-1,0] - u[:-2,0])) + \
	  0.5*Cy2*((q[1:-1,0] + q[1:-1,1])*(u[1:-1,1] - u[1:-1,0]) - \
	  (q[1:-1,0] + q[1:-1,1])*(u[1:-1,0] - u[1:-1,1])))
	  
    #y = Ly
    up[1:-1,-1] = c1*(2*u[1:-1,-1] - c2*um[1:-1,-1] + \
	  0.5*Cx2*((q[1:-1,-1] + q[2:,-1])*(u[2:,-1] - u[1:-1,-1]) - \
	  (q[1:-1,-1] + q[:-2,-1])*(u[1:-1,-1] - u[:-2,-1])) + \
	  0.5*Cy2*((q[1:-1,-1] + q[1:-1,-2])*(u[1:-1,-2] - u[1:-1,-1]) - \
	  (q[1:-1,-1] + q[1:-1,-2])*(u[1:-1,-1] - u[1:-1,-2])))
	  
    # Corners:
    up[0,0] = c1*(2*u[0,0] - c2*um[0,0] + \
            0.5*Cx2*((q[0,0] + q[1,0])*(u[1,0]-u[0,0]) - (q[0,0] + q[1,0])*(u[0,0] - u[1,0])) + \
            0.5*Cx2*((q[0,0] + q[0,1])*(u[0,1]-u[0,0]) - (q[0,0] + q[0,1])*(u[0,0] - u[0,1])))
    
    up[Nx,0] = c1*(2*u[Nx,0] - c2*um[Nx,0] + \
            0.5*Cx2*((q[Nx-1,0] + q[Nx,0])*(u[Nx-1,0] -u [Nx,0]) - (q[Nx,0] + q[Nx-1,0])*(u[Nx,0] - u[Nx-1,0])) + \
            0.5*Cx2*((q[Nx,0] + q[Nx,1])*(u[Nx,1]-u[Nx,0]) - (q[Nx,0] + q[Nx,1])*(u[Nx,0] - u[Nx,1])))
   
    up[0,Ny] = c1*(2*u[0,Ny] - c2*um[0,Ny] + \
            0.5*Cx2*((q[0,Ny] + q[1,Ny])*(u[1,Ny] - u [0,Ny]) - (q[1,Ny] + q[0,Ny])*(u[0,Ny] - u[1,Ny])) + \
            0.5*Cx2*((q[0,Ny] + q[0,Ny-1])*(u[0,Ny-1] - u[0,Ny]) - (q[0,Ny] + q[0,Ny-1])*(u[0,Ny] - u[0,Ny-1])))
   
    up[Nx,Ny] = c1*(2*u[Nx,Ny] - c2*um[Nx,Ny] + \
            0.5*Cx2*((q[Nx,Ny] + q[Nx-1,Ny])*(u[Nx-1,Ny] - u [Nx,Ny]) - (q[Nx-1,Ny] + q[Nx,Ny])*(u[Nx,Ny] - u[Nx-1,Ny])) + \
            0.5*Cx2*((q[Nx,Ny] + q[Nx,Ny-1])*(u[Nx,Ny-1] - u[Nx,Ny]) - (q[Nx,Ny] + q[Nx,Ny-1])*(u[Nx,Ny] - u[Nx,Ny-1])))
            
    
    um = u.copy()
    u = up.copy()
    return u,um

def plot_3D(x, y, u, s):
    #mlab.mesh(x,y,u)
    #mlab.surf(x,y,u)
    s.mlab_source.scalars = u[1:-1,1:-1]
    #mlab.savefig("wtmp%.5f.png" % n)

import glob, os
def run_Gaussian(version):
    dt = -1
    Nx = 50
    Ny = 50
    T = 10
    Lx = 10
    Ly = 10
    c = 1.1
    b = 0.1
    sigma = 0.8
    def I(x,y):
        """Gaussian peak at (Lx/2, Ly/2)."""
        return exp(-0.5*((x-Lx/2.0)/(sigma))**2 - 0.5*((y-Ly/2.0)/(sigma))**2)
    version = version
    V = ones((Nx+1,Ny+1))
    V = V*0.01
    q = ones((Nx+1,Ny+1))
    q = 0.5*q
    
    if version == "scalar":
        def f(x,y,t):
            return 0
    else:
        def f(x,y,t):
            return zeros((Nx+1,Ny+1))
    
    
    dt = solver(I, V, q, f, Lx, Ly, Nx, Ny, c, dt, T, b, version, make_plot=True)
    
    """
    movie("wtmp*.png")
    print "dt: ", dt
    for filename in glob("wtmp*.png"):
        os.remove(filename)
    """
    
    
run_Gaussian("vectorized")
#run_Gaussian("scalar")
print "No syntax errors. Yeey!"

   
