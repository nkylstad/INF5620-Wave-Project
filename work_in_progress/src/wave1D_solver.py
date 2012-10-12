from numpy import *
import math
import matplotlib.pyplot as plt
from mayavi import mlab

def solver_1D(Lx, Nx, T, dt, c, I, q, V, f, b, version):
    
    x = linspace(0,Lx,Nx+1)
    X = meshgrid(x,x)
    dx = float(x[1] - x[0]) 

    dt = c*dx
    N = int(round(float(T/dt)))
    t = linspace(0,T,N+1)
    c1 = 1/(1 + (b*dt)/2)
    c2 = 1 - (b*dt)/2
    Cx2 = (dt/dx)**2
    dt2 = dt**2
    
    u = zeros(Nx+1)  # The new soluion at the next timestep
    u_1 = zeros(Nx+1) # The solution from the current time step
    u_2 = zeros(Nx+1)  # The solution from the previous time step
    
    # Initial conditions
    if version == "scalar":
        for i in range(0,Nx+1):
            u_2[i] = I(x[i])
    else:  # vectorized version
        u_2[:] = I(X)
        
        
    # special scheme for the first step:
    if version == "scalar":
        for i in range(0,Nx+1):
            # Boundary conditions
            if i == 0: im1 = i+1
            else: im1 = i-1    # im1 represents the index i-1
            if i == Nx: ip1 = i-1
            else: ip1 = i+1  # ip1 represents the index i+1
                
            u_1[i] = u_2[i] + dt*c2*V(x[i]) + \
                0.5*Cx2*((q[ip1,j] + q[i])*(u_2[ip1] - u_2[i]) - (q[i] + q[im1])*(u_2[i] - u_2[im1])) + \
                dt2*f(x[i],t[0])
                            
    else:  #vectorized version
        u_1[1:-1] = u_2[1:-1] + dt*c2*V[1:-1] + \
            0.5*Cx2*((q[2:] + q[1:-1])*(u_2[2:] - u_2[1:-1]) - (q[1:-1] + q[:-2])*(u_2[1:-1] - u_2[:-2])) + \
            
        # boundary conditions:
        u_1[0,:] = u_1[1,:]
        u_1[-1,:] = u_1[-2,:]
   
    E = zeros(N+1)

    for n in range(2,N+1):
        
        if version == "scalar":
            u_1,u_2 = advance_scalar(u, u_1, u_2, Nx, x, q, f, c1, c2, Cx2, dt2, b, t[n])
                plot_1D(u_1,x)  
        else:
            u_1,u_2 = advance_vectorized(u, u_1, u_2, q, f, c1, c2, Cx2, Cy2, t)
            #plot_u(u_1, x, X, y, Y, t, n, 2)
       
    return u_1       
         
         
def advance_scalar_1D(u, u_1, u_2, Nx, x, q, f, c1, c2, Cx2, dt2, b, tn):
    for i in range(0,Nx+1):
        # Boundary conditions
        if i == 0: im1 = i+1
        else: im1 = i-1    # im1 represents the index i-1
        if i == Nx: ip1 = i-1
        else: ip1 = i+1  # ip1 represents the index i+1
        if j == 0: jm1 = j+1
        else: jm1 = j-1    # jm1 represents the index j-1
        if j == Ny: jp1 = j-1
        else: jp1 = j+1  # jp1 represents the index j+1
        u[i] = c1*(2*u_1[i] - c2*u_2[i] + \
            0.5*Cx2*((q[ip1] + q[i])*(u_2[ip1] - u_2[i]) - (q[i] + q[im1])*(u_2[i] - u_2[im1])) + \
            dt2*f(x[i],tn))
                
    u_2 = u_1.copy()
    u_1 = u.copy()
    return u_1, u_2 
    