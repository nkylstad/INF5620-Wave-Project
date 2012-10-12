from numpy import *
import math, os, subprocess, sys, glob
from plot_u import plot_u, make_movie

def solver(Lx, Ly, Nx, Ny, T, dt, c, I, q, V, f, b, version, h=0.4, oneD=False, standing=False, make_plot=True):
    
    x = linspace(0,Lx,Nx+1)
    y = linspace(0,Ly,Ny+1)
    X,Y = meshgrid(x,y)
    dx = float(x[1] - x[0]) 
    dy = float(y[1] - y[0])
   
    stability_limit = (1/float(c))*(1/sqrt(1/dx**2 + 1/dy**2))
    if dt <= 0:
        dt = 1*stability_limit
    elif dt > stability_limit:
        print "Error: dt too large."
    
    if oneD:
        y = zeros(Ny+1)
        dy = 1
        dt = dx/c
        Y = meshgrid(y, y)
    
    if standing:
        Fx = 0.8
        Fy = 0.8
        Ft = 1/c *1/sqrt(1/Fx**2 + 1/Fy**2)
        dx = Fx*h
        dy = Fy*h
        dt = Ft*h
    
    N = int(round(float(T/dt)))
    t = linspace(0,T,N+1)
    tv = meshgrid(t,t)
    c1 = 1/(1 + (b*dt)/2)
    c2 = 1 - (b*dt)/2
    Cx2 = (dt/dx)**2
    Cy2 = (dt/dy)**2
    dt2 = dt**2
    
    if oneD:
        Cy2 = 0
        
    file = open("initial.txt",'w')
    file.write("Nx="+str(Nx)+"\n")
    file.write("Lx="+str(Lx)+"\n")
    file.write("Ny="+str(Ny)+"\n")
    file.write("Ly="+str(Ly)+"\n")
    file.write("T="+str(T)+"\n")
    file.write("N="+str(N)+"\n")
    file.close()
    
    u = zeros((Nx+1, Ny+1))  # The new soluion at the next timestep
    u_1 = zeros((Nx+1, Ny+1)) # The solution from the current time step
    u_2 = zeros((Nx+1, Ny+1))  # The solution from the previous time step
    
    # Initial conditions
    if version == "scalar":
        for i in range(0,Nx+1):
            for j in range(0,Ny+1):
                u_2[i,j] = I(x[i],y[j])
    else:  # vectorized version
        u_2[:,:] = I(X,Y)
        

    
    if make_plot:
        if oneD:
            plot_u(u_2, x, t[0], 0, b, Lx)
        else:
            savetxt("u0.txt", u_2)
        
    # special scheme for the first step:
    if version == "scalar":
        for i in range(0,Nx+1):
            for j in range(0,Ny+1):
                # Boundary conditions
                if i == 0: im1 = i+1
                else: im1 = i-1    # im1 represents the index i-1
                if i == Nx: ip1 = i-1
                else: ip1 = i+1  # ip1 represents the index i+1
                if j == 0: jm1 = j+1
                else: jm1 = j-1    # jm1 represents the index j-1
                if j == Ny: jp1 = j-1
                else: jp1 = j+1  # jp1 represents the index j+1
                
                u_1[i,j] = u_2[i,j] + dt*c2*V(x[i],y[j]) + \
                            0.25*Cx2*((q[ip1,j] + q[i,j])*(u_2[ip1,j] - u_2[i,j]) - (q[i,j] + q[im1,j])*(u_2[i,j] - u_2[im1,j])) + \
                            0.25*Cy2*((q[i,jp1] + q[i,j])*(u_2[i,jp1] - u_2[i,j]) - (q[i,j] + q[i,jm1])*(u_2[i,j] - u_2[i,jm1])) + \
                            0.5*dt2*f(x[i],y[j],t[0])
                            
    else:  #vectorized version
        Vv = V(X,Y)
        fv = f(X,Y,t[0])
        u_1[1:-1,1:-1] = u_2[1:-1,1:-1] + dt*c2*Vv[1:-1,1:-1] + \
            0.5**2*Cx2*((q[2:,1:-1] + q[1:-1,1:-1])*(u_2[2:,1:-1] - u_2[1:-1,1:-1]) - (q[1:-1,1:-1] + q[:-2,1:-1])*(u_2[1:-1,1:-1] - u_2[:-2,1:-1])) + \
            0.5**2*Cy2*((q[1:-1,2:] + q[1:-1,1:-1])*(u_2[1:-1,2:] - u_2[1:-1,1:-1]) - (q[1:-1,1:-1] + q[1:-1,:-2])*(u_2[1:-1,1:-1] - u_2[1:-1,:-2])) \
            + 0.5*dt2*fv[1:-1,1:-1]
            
        # boundary conditions:
        u_1[0,:] = u_1[1,:]
        u_1[-1,:] = u_1[-2,:]
        u_1[:,0] = u_1[:,1]
        u_1[:,-1] = u_1[:,-2]
     
    if make_plot:
        if oneD:
            plot_u(u_1, x, t[1], 1, b, Lx)
        else:
            savetxt("texttmp%.4d.txt"%1, u_1)
    E = 0
    E_val = -1

    for n in range(1,N+1):
        if version == "scalar":
            u_1,u_2, E_val = advance_scalar(u, u_1, u_2, Nx, Ny, x, y, q, f, c1, c2, Cx2, Cy2, dt2, b, h, t[n], exact_manufactured_solution=None)
            if make_plot:
                if oneD:
                    plot_u(u_1, x, t[n], n, b, Lx)
                else:
                    if n%3 == 0:
                        savetxt("texttmp%.4d.txt" %n, u_1)
                    
        else:
            fv = f(X,Y,t[n])
            u_1,u_2 = advance_vectorized(u, u_1, u_2, q, fv, c1, c2, Cx2, Cy2, dt2, t, exact_manufactured_solution=None)
            if make_plot:
                if oneD:
                    plot_u(u_1, x, t[n], n, b, Lx)
                else:
                    if n%3 == 0:
                        savetxt("texttmp%.4d.txt" %n, u_1)
        
        if E < E_val:
            E = E_val
        
    
    #print "dt", dt 
    #print "E: ", E
    #print "E=dt^2+dx^2+dy^2:", dt**2+dx**2+dy**2
    return E     
         
         
def advance_scalar(u, u_1, u_2, Nx, Ny, x, y, q, f, c1, c2, Cx2, Cy2, dt2, b, h, tn, exact):
    error = 0
    for i in range(0,Nx+1):
        for j in range(0,Ny+1):
            
            # Boundary conditions
            if i == 0: im1 = i+1
            else: im1 = i-1    # im1 represents the index i-1
            if i == Nx: ip1 = i-1
            else: ip1 = i+1  # ip1 represents the index i+1
            if j == 0: jm1 = j+1
            else: jm1 = j-1    # jm1 represents the index j-1
            if j == Ny: jp1 = j-1
            else: jp1 = j+1  # jp1 represents the index j+1
            u[i,j] = c1*(2*u_1[i,j] - c2*u_2[i,j] + \
                0.5*Cx2*((q[ip1,j] + q[i,j])*(u_1[ip1,j] - u_1[i,j]) - (q[i,j] + q[im1,j])*(u_1[i,j] - u_1[im1,j])) + \
                0.5*Cy2*((q[i,jp1] + q[i,j])*(u_1[i,jp1] - u_1[i,j]) - (q[i,j] + q[i,jm1])*(u_1[i,j] - u_1[i,jm1])) + \
                dt2*f(x[i],y[j],tn))
            if exact:    
				u_exact = exact(x[i],y[j],b,tn)
				error = u_exact - u
            
    E_val = sqrt(dt*sum(error**2))
    #print "E_val:", E_val
                
    u_2 = u_1.copy()
    u_1 = u.copy()
    return u_1, u_2, E_val
    
def advance_vectorized(u, u_1, u_2, q, f, c1, c2, Cx2, Cy2, dt2, t, exact_manufactured_solution):
    u[1:-1,1:-1] = c1*(2*u_1[1:-1,1:-1] - c2*u_2[1:-1,1:-1] + \
      0.5*Cx2*((q[1:-1,1:-1] + q[2:,1:-1])*(u_1[2:,1:-1]-u_1[1:-1,1:-1]) - \
      (q[1:-1,1:-1] + q[:-2,1:-1])*(u_1[1:-1,1:-1] - u_1[:-2,1:-1])) + \
      0.5*Cy2*((q[1:-1,1:-1] + q[1:-1,2:])*(u_1[1:-1,2:] - u_1[1:-1,1:-1]) - \
      (q[1:-1,1:-1] + q[1:-1,:-2])*(u_1[1:-1,1:-1] - u_1[1:-1,:-2])) + dt2*f[1:-1,1:-1])
      
    #boundary conditions:
    
    u[0,:] = u[1,:]
    u[-1,:] = u[-2,:]
    u[:,0] = u[:,1]
    u[:,-1] = u[:,-2]
    
    u_2 = u_1.copy()
    u_1 = u.copy()
    return u_1, u_2





    
def test_1D_plug(version):
    C = 1
    Nx = 50
    Ny = 50
    Lx = 1
    Ly = 1
    T = 10
    c = 1
    b = 0.1
    dt = -1.0
    sigma = Lx/10.
    xc = Lx/2.
    #xc = 0
    
    q = ones((Nx+1,Ny+1))
    
    if version == "scalar":
        def V(x,y):
            return
    else: 
        V = zeros((Nx+1, Ny+1))
    
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
        
    if version == "scalar":
        def f(x,y,t):
            return 0
    else:
        f = zeros((Nx+1, Ny+1))
        
    u = solver(Lx, Ly, Nx, Ny, T, dt, c, I, q, V, f, b, version, oneD=True)
    

#test_constant_1D("scalar")    
#test_1D_plug("vectorized")



    
