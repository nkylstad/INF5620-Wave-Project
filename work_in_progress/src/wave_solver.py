from numpy import *
import math, os, subprocess, sys, glob
from plot_u import plot_u, make_movie

def solver(Lx, Ly, Nx, Ny, T, dt, c, I, q, V, f, b, version, h=0.4, w=None, exact=None, oneD=False, standing=False, make_plot=True):
    
    x = linspace(0,Lx,Nx+1)
    y = linspace(0,Ly,Ny+1)
    X,Y = meshgrid(x,y)     # Create spatial points
    dx = float(x[1] - x[0])   # Calculate dx
    dy = float(y[1] - y[0])   # Calculate dy
   
    stability_limit = (1/float(c))*(1/sqrt(1/dx**2 + 1/dy**2))  # optimal value for dt
    if dt <= 0:
        dt = 1*stability_limit
    elif dt > stability_limit:
        print "Error: dt too large."
    
    if oneD:    # Special case for 1D
        y = zeros(Ny+1)
        dy = 1
        dt = dx/c
        Y = meshgrid(y, y)
    
    if standing:  # Special case for standing wave (to be used in verification test)
        Fx = 0.8
        Fy = 0.8
        Ft = 1/c *1/sqrt(1/Fx**2 + 1/Fy**2)
        dx = Fx*h
        dy = Fy*h
        dt = Ft*h
        Nx = int(round(float(Lx/dx)))
        Ny = int(round(float(Ly/dy)))
        x = linspace(0, Lx, Nx+1)
        y = linspace(0, Ly, Ny+1)
        X, Y = meshgrid(x,y)
    
    N = int(round(float(T/dt)))
    t = linspace(0,T,N+1)
    c1 = 1/(1 + (b*dt)/2)  # Set constants / help variables to be used in the scheme:
    c2 = 1 - (b*dt)/2
    Cx2 = (dt/dx)**2
    Cy2 = (dt/dy)**2
    dt2 = dt**2
    
    if oneD:
        Cy2 = 0  # Special case for 1D
        
    if make_plot:    
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
    q = q(X,Y)
    E_list = zeros(N)
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
                
                # Scheme for all points (including boundary)
                u_1[i,j] = u_2[i,j] + dt*c2*V(x[i],y[j]) + \
                            0.25*Cx2*((q[ip1,j] + q[i,j])*(u_2[ip1,j] - u_2[i,j]) - (q[i,j] + q[im1,j])*(u_2[i,j] - u_2[im1,j])) + \
                            0.25*Cy2*((q[i,jp1] + q[i,j])*(u_2[i,jp1] - u_2[i,j]) - (q[i,j] + q[i,jm1])*(u_2[i,j] - u_2[i,jm1])) + \
                            0.5*dt2*f(x[i],y[j],t[0])
                                                    
    else:  #vectorized version
        Vv = V(X,Y)
        fv = f(X,Y,t[0])
        #print fv
        # Scheme for all interior points
        u_1[1:-1,1:-1] = u_2[1:-1,1:-1] + dt*c2*Vv[1:-1,1:-1] + \
            0.5**2*Cx2*((q[2:,1:-1] + q[1:-1,1:-1])*(u_2[2:,1:-1] - u_2[1:-1,1:-1]) - (q[1:-1,1:-1] + q[:-2,1:-1])*(u_2[1:-1,1:-1] - u_2[:-2,1:-1])) + \
            0.5**2*Cy2*((q[1:-1,2:] + q[1:-1,1:-1])*(u_2[1:-1,2:] - u_2[1:-1,1:-1]) - (q[1:-1,1:-1] + q[1:-1,:-2])*(u_2[1:-1,1:-1] - u_2[1:-1,:-2])) \
            + 0.5*dt2*fv[1:-1,1:-1]
            
        # simple boundary conditions:
        u_1[0,:] = u_1[1,:]
        u_1[-1,:] = u_1[-2,:]
        u_1[:,0] = u_1[:,1]
        u_1[:,-1] = u_1[:,-2]
        
        # boundary conditions
        #u_1[0,1:-1] = u_2[0,1:-1] + dt*c2*Vv[0,1:-1] + \
                #0.5**2*Cx2*((q[1,1:-1] + q[0,1:-1])*(u_2[1,1:-1] - u_2[0,1:-1]) - (q[0,1:-1] + q[1,1:-1])*(u[0,1:-1] - u[1,1:-1])) + \
                #0.5**2*Cy2*((q[0,2:] + q[0,1:-1])*(u_2[0,2:] - u_2[0,1:-1]) - (q[0,1:-1] + q[0,:-2])*(u[0,1:-1] - u[0,:-2])) + \
                #0.5*dt2*fv[0,1:-1]
                
        #u_1[-1,1:-1] = u_2[-1,1:-1] + dt*c2*Vv[-1,1:-1] + \
                #0.5**2*Cx2*((q[-2,1:-1] + q[-1,1:-1])*(u_2[-2,1:-1] - u_2[-1,1:-1]) - (q[-1,1:-1] + q[-2,1:-1])*(u[-1,1:-1] - u[-2,1:-1])) + \
                #0.5**2*Cy2*((q[-1,2:] + q[-1,1:-1])*(u_2[-1,2:] - u_2[-1,1:-1]) - (q[-1,1:-1] + q[-1,:-2])*(u[-1,1:-1] - u[-1,:-2])) + \
                #0.5*dt2*fv[-1,1:-1]
                
        #u_1[1:-1,0] = u_2[1:-1, 0] + dt*c2*Vv[1:-1, 0] + \
                #0.5**2*Cx2*((q[2:,0] + q[1:-1,0])*(u_2[2:,0] - u_2[1:-1,0]) - (q[1:-1,0] + q[:-2,0])*(u[1:-1,0] - u[:-2,0])) + \
                #0.5**2*Cy2*((q[1:-1,1] + q[1:-1,0])*(u_2[1:-1,1] - u_2[1:-1,0]) - (q[1:-1,0] + q[1:-1,1])*(u[1:-1,0] - u[1:-1,1])) + \
                #0.5*dt2*fv[1:-1,0]
                
        #u_1[1:-1,-1] = u_2[1:-1,-1] + dt*c2*Vv[1:-1,-1] + \
                #0.5**2*Cx2*((q[2:,-1] + q[1:-1,-1])*(u_2[2:,-1] - u_2[1:-1,-1]) - (q[1:-1,-1] + q[:-2,-1])*(u[1:-1,-1] - u[:-2,-1])) + \
                #0.5**2*Cy2*((q[1:-1,-2] + q[1:-1,-1])*(u_2[1:-1,-2] - u_2[1:-1,-1]) - (q[1:-1,-1] + q[1:-1,-2])*(u[1:-1,-1] - u[1:-1,-2])) + \
                #0.5*dt2*fv[1:-1,-1]
                
        ## Corners:
        #u_1[0,0] = u_2[0,0] + dt*c2*Vv[0,0] + \
                #0.5**2*Cx2*((q[1,0] + q[0,0])*(u_2[1,0] - u_2[0,0]) - (q[0,0] + q[1,0])*(u[0,0] - u[1,0])) + \
                #0.5**2*Cy2*((q[0,1] + q[0,0])*(u_2[0,1] - u_2[0,0]) - (q[0,0] + q[0,1])*(u[0,0] - u[0,1])) + \
                #0.5*dt2*fv[0,0]
                
        #u_1[0,-1] = u_2[0,-1] + dt*c2*Vv[0,-1] + \
                #0.5**2*Cx2*((q[1,-1] + q[0,-1])*(u_2[1,-1] - u_2[0,-1]) - (q[0,-1] + q[1,-1])*(u[0,-1] - u[1,-1])) + \
                #0.5**2*Cy2*((q[0,-2] + q[0,-1])*(u_2[0,-2] - u_2[0,-1]) - (q[0,-1] + q[0,-2])*(u[0,-1] - u[0,-2])) + \
                #0.5*dt2*fv[0,-1]
                
        #u_1[-1,0] = u_2[-1,0] + dt*c2*Vv[-1,0] + \
                #0.5**2*Cx2*((q[-2,0] + q[-1,0])*(u_2[-2,0] - u_2[-1,0]) - (q[-1,0] + q[-2,0])*(u[-1,0] - u[-2,0])) + \
                #0.5**2*Cy2*((q[-1,1] + q[-1,0])*(u_2[-1,1] - u_2[-1,0]) - (q[-1,0] + q[-1,1])*(u[-1,0] - u[-1,1])) + \
                #0.5*dt2*fv[0,0]
                
        #u_1[-1,-1] = u_2[-1,-1] + dt*c2*Vv[-1,-1] + \
                #0.5**2*Cx2*((q[-2,-1] + q[-1,-1])*(u_2[-2,-1] - u_2[-1,-1]) - (q[-1,-1] + q[-2,-1])*(u[-1,-1] - u[-2,-1])) + \
                #0.5**2*Cy2*((q[-1,-2] + q[-1,-1])*(u_2[-1,-2] - u_2[-1,-1]) - (q[-1,-1] + q[-1,-2])*(u[-1,-1] - u[-1,-2])) + \
                #0.5*dt2*fv[0,-1]
                
    if exact:    
        Err = exact(X, Y, b, t[1], w) - u_1   # compare numerical solution to exact solution
        E = sqrt(h*sum(Err**2))
        E_list[0] = E
    
    if make_plot:
        if oneD:
            plot_u(u_1, x, t[1], 1, b, Lx)
            make_movie()
        else:
            savetxt("texttmp%.4d.txt"%1, u_1)
   

    for n in range(1,N):
        if version == "scalar":
            u_1,u_2,E = advance_scalar(u, u_1, u_2, Nx, Ny, x, y, q, f, c1, c2, Cx2, Cy2, dt2, b, h, w, t[n+1], exact)
            if make_plot:
                if oneD:
                    plot_u(u_1, x, t[n], n, b, Lx)
                    make_movie()
                else:
                    #if n%3 == 0:
                    savetxt("texttmp%.4d.txt" %n, u_1)
                    
        else:   # vectorized version
            fv = f(X,Y,t[n])
            u_1,u_2,E = advance_vectorized(X, Y, u, u_1, u_2, q, fv, c1, c2, Cx2, Cy2, dt2, t[n+1], b, h, w, Nx, Ny, exact)
            if make_plot:
                if oneD:
                    plot_u(u_1, x, t[n], n, b, Lx)
                else:
                    #if n%3 == 0:
                    savetxt("texttmp%.4d.txt" %n, u_1)
        
        if exact:
            E_list[n] = E
            
    return E_list, u_1
         
         
         
def advance_scalar(u, u_1, u_2, Nx, Ny, x, y, q, f, c1, c2, Cx2, Cy2, dt2, b, h, w, tn, exact):
    Err = zeros((Nx+1,Ny+1))
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
            
            # Scheme for all points (including boundary)
            u[i,j] = c1*(2*u_1[i,j] - c2*u_2[i,j] + \
                0.5*Cx2*((q[ip1,j] + q[i,j])*(u_1[ip1,j] - u_1[i,j]) - (q[i,j] + q[im1,j])*(u_1[i,j] - u_1[im1,j])) + \
                0.5*Cy2*((q[i,jp1] + q[i,j])*(u_1[i,jp1] - u_1[i,j]) - (q[i,j] + q[i,jm1])*(u_1[i,j] - u_1[i,jm1])) + \
                dt2*f(x[i],y[j],tn))
            
            if exact:
                Err[i,j] = exact(x[i],y[i],b,tn,w,version="scalar") - u[i,j]  # compare numerical solution with exact solution
    
    E = sqrt(h*sum(Err**2))
    u_2 = u_1.copy()
    u_1 = u.copy()
    return u_1, u_2, E
    
    
    
def advance_vectorized(X, Y ,u, u_1, u_2, q, f, c1, c2, Cx2, Cy2, dt2, tn, b, h, w, Nx, Ny, exact):
    
    # Scheme for all interior points
    u[1:-1,1:-1] = c1*(2*u_1[1:-1,1:-1] - c2*u_2[1:-1,1:-1] + \
      0.5*Cx2*((q[1:-1,1:-1] + q[2:,1:-1])*(u_1[2:,1:-1]-u_1[1:-1,1:-1]) - \
      (q[1:-1,1:-1] + q[:-2,1:-1])*(u_1[1:-1,1:-1] - u_1[:-2,1:-1])) + \
      0.5*Cy2*((q[1:-1,1:-1] + q[1:-1,2:])*(u_1[1:-1,2:] - u_1[1:-1,1:-1]) - \
      (q[1:-1,1:-1] + q[1:-1,:-2])*(u_1[1:-1,1:-1] - u_1[1:-1,:-2])) + dt2*f[1:-1,1:-1])
      
    # simple boundary conditions:
    u[0,:] = u[1,:]
    u[-1,:] = u[-2,:]
    u[:,0] = u[:,1]
    u[:,-1] = u[:,-2]
    
    # boundary conditions:
    # x = 0
    #u[0,1:-1] = c1*(2*u_1[0,1:-1] - c2*u_2[0,1:-1] + \
      #0.5*Cx2*((q[0,1:-1] + q[1,1:-1])*(u_1[1,1:-1] - u_1[0,1:-1]) - \
      #(q[0,1:-1] + q[1,1:-1])*(u_1[0,1:-1] - u_1[1,1:-1])) + \
      #0.5*Cy2*((q[0,1:-1] + q[0,2:])*(u_1[0,2:] - u_1[0,1:-1]) - \
      #(q[0,1:-1] + q[0,:-2])*(u_1[0,1:-1] - u_1[0,:-2])) + dt2*f[0, 1:-1])
    
    ## x = Lx
    #u[-1,1:-1] = c1*(2*u_1[-1,1:-1] - c2*u_2[-1,1:-1] + \
      #0.5*Cx2*((q[-1,1:-1] + q[-2,1:-1])*(u_1[-2,1:-1] - u_1[-1,1:-1]) - \
      #(q[-1,1:-1] + q[-2,1:-1])*(u_1[-1,1:-1] - u_1[-2,1:-1])) + \
      #0.5*Cy2*((q[-1,1:-1] + q[-1,2:])*(u_1[-1,2:] - u_1[-1,1:-1]) - \
      #(q[-1,1:-1] + q[-1,:-2])*(u_1[-1,1:-1] - u_1[-1,:-2])) + dt2*f[-1, 1:-1])

    ##y = 0
    #u[1:-1,0] = c1*(2*u_1[1:-1,0] - c2*u_2[1:-1,0] + \
      #0.5*Cx2*((q[1:-1,0] + q[2:,0])*(u_1[2:,0] - u_1[1:-1,0]) - \
      #(q[1:-1,0] + q[:-2,0])*(u_1[1:-1,0] - u_1[:-2,0])) + \
      #0.5*Cy2*((q[1:-1,0] + q[1:-1,1])*(u_1[1:-1,1] - u_1[1:-1,0]) - \
      #(q[1:-1,0] + q[1:-1,1])*(u_1[1:-1,0] - u_1[1:-1,1])) + dt2*f[1:-1, 0])
      
    ##y = Ly
    #u[1:-1,-1] = c1*(2*u_1[1:-1,-1] - c2*u_2[1:-1,-1] + \
      #0.5*Cx2*((q[1:-1,-1] + q[2:,-1])*(u_1[2:,-1] - u_1[1:-1,-1]) - \
      #(q[1:-1,-1] + q[:-2,-1])*(u_1[1:-1,-1] - u_1[:-2,-1])) + \
      #0.5*Cy2*((q[1:-1,-1] + q[1:-1,-2])*(u_1[1:-1,-2] - u_1[1:-1,-1]) - \
      #(q[1:-1,-1] + q[1:-1,-2])*(u_1[1:-1,-1] - u_1[1:-1,-2])) + dt2*f[1:-1,-1])
      
    ## Corners:
    ##x=0, y=0
    #u[0,0] = c1*(2*u_1[0,0] - c2*u_2[0,0] + \
            #0.5*Cx2*((q[0,0] + q[1,0])*(u_1[1,0]-u_1[0,0]) - (q[0,0] + q[1,0])*(u_1[0,0] - u_1[1,0])) + \
            #0.5*Cx2*((q[0,0] + q[0,1])*(u_1[0,1]-u_1[0,0]) - (q[0,0] + q[0,1])*(u_1[0,0] - u_1[0,1])) + dt2*f[0, 0])
    
    ## x=Lx, y=0
    #u[Nx,0] = c1*(2*u_1[Nx,0] - c2*u_2[Nx,0] + \
            #0.5*Cx2*((q[Nx-1,0] + q[Nx,0])*(u_1[Nx-1,0] -u_1 [Nx,0]) - (q[Nx,0] + q[Nx-1,0])*(u_1[Nx,0] - u_1[Nx-1,0])) + \
            #0.5*Cx2*((q[Nx,0] + q[Nx,1])*(u_1[Nx,1]-u_1[Nx,0]) - (q[Nx,0] + q[Nx,1])*(u_1[Nx,0] - u_1[Nx,1])) + dt2*f[Nx, 0])
  
  ## x=0, y=Ly
    #u[0,Ny] = c1*(2*u_1[0,Ny] - c2*u_2[0,Ny] + \
            #0.5*Cx2*((q[0,Ny] + q[1,Ny])*(u_1[1,Ny] - u_1 [0,Ny]) - (q[1,Ny] + q[0,Ny])*(u_1[0,Ny] - u_1[1,Ny])) + \
            #0.5*Cx2*((q[0,Ny] + q[0,Ny-1])*(u_1[0,Ny-1] - u_1[0,Ny]) - (q[0,Ny] + q[0,Ny-1])*(u_1[0,Ny] - u_1[0,Ny-1])) + dt2*f[0, Ny])
   
   ## x=Lx, y=Ly
    #u[Nx,Ny] = c1*(2*u_1[Nx,Ny] - c2*u_2[Nx,Ny] + \
            #0.5*Cx2*((q[Nx,Ny] + q[Nx-1,Ny])*(u_1[Nx-1,Ny] - u_1 [Nx,Ny]) - (q[Nx-1,Ny] + q[Nx,Ny])*(u_1[Nx,Ny] - u_1[Nx-1,Ny])) + \
            #0.5*Cx2*((q[Nx,Ny] + q[Nx,Ny-1])*(u_1[Nx,Ny-1] - u_1[Nx,Ny]) - (q[Nx,Ny] + q[Nx,Ny-1])*(u_1[Nx,Ny] - u_1[Nx,Ny-1])) + dt2*f[Nx, Ny])
    
    if exact:
        Err = exact(X,Y,b,tn,w) - u  # Compare numerical soluion with exact solution
        E = sqrt(h*sum(Err**2))
    else:
        E = 0
    
    u_2 = u_1.copy()
    u_1 = u.copy()
    return u_1, u_2, E
