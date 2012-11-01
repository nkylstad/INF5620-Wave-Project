from numpy import *
import math, os, subprocess, sys, glob
from plot_u import plot_u, make_movie

def solver(Lx, Ly, Nx, Ny, T, dt, c, I, q, V, f, b, version, B=None, w=None, exact=None, oneD=False, make_plot=True):
    
    dx_x = linspace(0, Lx, Nx+1)
    dy_y = linspace(0, Ly, Ny+1)

    if version == "scalar":
        x = linspace(0, Lx, Nx+1)
        y = linspace(0, Ly, Ny+1)
    else:
        x = linspace(0, Lx, Nx+3)
        y = linspace(0, Ly, Ny+3)
    X,Y = meshgrid(x,y)     # Create spatial points
    dx = float(dx_x[1] - dx_x[0])   # Calculate dx
    dy = float(dy_y[1] - dy_y[0])   # Calculate dy

    q = q(X,Y)
    q_max = amax(q)

    if make_plot:
        if B:
            B = B(X,Y)
            if version == "scalar":
                savetxt("hill.txt", B)
            else:
                savetxt("hill.txt", B[1:-1,1:-1])

   
    stability_limit = (1/float(sqrt(q_max)))*(1/sqrt(1/dx**2 + 1/dy**2))  # optimal value for dt
    if dt <= 0:
        dt = 1*stability_limit
    elif dt > stability_limit:
        print "Error: dt too large."

    if oneD:    # Special case for 1D
        y = y*0
        dy = 1
        dt = dx/c
        Y,Y = meshgrid(y,y)
        print "bleh"
        
    print "dx: ", dx
    print "dy: ", dy

    
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
    
    if version == "scalar":
        u = zeros((Nx+1, Ny+1))  # The new soluion at the next timestep
        u_1 = zeros((Nx+1, Ny+1)) # The solution from the current time step
        u_2 = zeros((Nx+1, Ny+1))  # The solution from the previous time step
    else:
        u = zeros((Nx+3, Ny+3))  # The new soluion at the next timestep
        u_1 = zeros((Nx+3, Ny+3)) # The solution from the current time step
        u_2 = zeros((Nx+3, Ny+3))  # The solution from the previous time step
    
    # Initial conditions
    if version == "scalar":
        for i in range(0,Nx+1):
            for j in range(0,Ny+1):
                u_2[i,j] = I(x[i],y[j])
    else:  # vectorized version
        u_2[:,:] = I(X,Y)
    
    
    if make_plot:
        if oneD:
            if version == "scalar":
                plot_u(u_2, x, t[0], 0, b, Lx)
            else:
                plot_u(u_2[1:-1,1:-1], x[1:-1], t[0], 0, b, Lx)
        else:
            if version == "scalar":
                savetxt("u0.txt", u_2)
            else:
                savetxt("u0.txt", u_2[1:-1,1:-1])
   
    Vv = V(X,Y)
    fv = f(X,Y,t[0])
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
                u_1[i,j] = u_2[i,j] + dt*c2*Vv[i,j] + \
                            0.5**2*Cx2*((q[ip1,j] + q[i,j])*(u_2[ip1,j] - u_2[i,j]) - (q[i,j] + q[im1,j])*(u_2[i,j] - u_2[im1,j])) + \
                            0.5**2*Cy2*((q[i,jp1] + q[i,j])*(u_2[i,jp1] - u_2[i,j]) - (q[i,j] + q[i,jm1])*(u_2[i,j] - u_2[i,jm1])) + \
                            0.5*dt2*fv[i,j]
                                                    
    else:  #vectorized version
        
        # boundary conditions:
        u_2[0,:] = u_2[2,:] # u[-1, y] = u[1, y]
        u_2[-1,:] = u_2[-3,:] # u[Nx+1, y] = u[Nx-1, y]
        u_2[:,0] = u_2[:,2] # u[x, -1] = u[x, 1]
        u_2[:,-1] = u_2[:,-3] # u[x, Ny+1] = u[x, Ny-1]
        
        # Scheme for all interior points
        u_1[1:-1,1:-1] = u_2[1:-1,1:-1] + dt*c2*Vv[1:-1,1:-1] + \
            0.5**2*Cx2*((q[2:,1:-1] + q[1:-1,1:-1])*(u_2[2:,1:-1] - u_2[1:-1,1:-1]) - (q[1:-1,1:-1] + q[:-2,1:-1])*(u_2[1:-1,1:-1] - u_2[:-2,1:-1])) + \
            0.5**2*Cy2*((q[1:-1,2:] + q[1:-1,1:-1])*(u_2[1:-1,2:] - u_2[1:-1,1:-1]) - (q[1:-1,1:-1] + q[1:-1,:-2])*(u_2[1:-1,1:-1] - u_2[1:-1,:-2])) \
            + 0.5*dt2*fv[1:-1,1:-1]
                
    if exact is not None:  
        exact_v = exact(X,Y,b,t[1],w)  
        if version=="scalar":
            Err = exact_v - u_1   # compare numerical solution to exact solution
        else:
            Err = exact_v[1:-1,1:-1] - u_1[1:-1,1:-1]
        E = sqrt(dx*dx*sum(Err**2))
        E_list[0] = E
    else:
        exact_v = None
    
    if make_plot:
        if oneD:
            if version == "scalar":
                plot_u(u_1, x, t[1], 1, b, Lx)
            else:
                plot_u(u_1[1:-1,1:-1], x[1:-1], t[1], 1, b, Lx)
        else:
            if version == "scalar":
                savetxt("texttmp%.4d.txt"%1, u_1)
            else:
                savetxt("texttmp%.4d.txt"%1, u_1[1:-1, 1:-1])
   

    for n in range(1,N):
        fv = f(X,Y,t[n])
        if exact is not None:
            exact_v = exact(X,Y,b,t[n+1],w)
        if version == "scalar":
            u_1, u_2, E = advance_scalar(u, u_1, u_2, Nx, Ny, x, y, q, fv, c1, c2, Cx2, Cy2, dt2, b, dx, w, t[n+1], exact_v)
            if make_plot:
                if oneD:
                    plot_u(u_1, x, t[n+1], n, b, Lx)
                else:
                    #if n%3 == 0:
                    savetxt("texttmp%.4d.txt" %n, u_1)
                    
        else:   # vectorized version
            u_1, u_2, E = advance_vectorized(X, Y, u, u_1, u_2, q, fv, c1, c2, Cx2, Cy2, dt2, t[n+1], b, dx, w, Nx, Ny, exact_v)
            if make_plot:
                if oneD:
                    plot_u(u_1[1:-1,1:-1], x[1:-1], t[n+1], n, b, Lx)
                else:
                    #if n%3 == 0:
                    savetxt("texttmp%.4d.txt" %n, u_1[1:-1, 1:-1])
        
        if exact:
            E_list[n] = E
    if make_plot:
        if oneD:
            make_movie()        
    return E_list, u_1, dx
         
         
         
def advance_scalar(u, u_1, u_2, Nx, Ny, x, y, q, f, c1, c2, Cx2, Cy2, dt2, b, dx, w, tn, exact):
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
                dt2*f[i,j])
            
    if exact is not None:
        Err = exact - u  # compare numerical solution with exact solution
        E = sqrt(dx*dx*sum(Err**2))
    else:
        E = 0

    u_2 = u_1.copy()
    u_1 = u.copy()
    return u_1, u_2, E
    
    
    
def advance_vectorized(X, Y ,u, u_1, u_2, q, f, c1, c2, Cx2, Cy2, dt2, tn, b, dx, w, Nx, Ny, exact):
    # boundary conditions:
    u_1[0,:] = u_1[2,:] # u[-1, y] = u[1, y]
    u_1[-1,:] = u_1[-3,:] # u[Nx+1, y] = u[Nx-1, y]
    u_1[:,0] = u_1[:,2] # u[x, -1] = u[x, 1]
    u_1[:,-1] = u_1[:,-3] # u[x, Ny+1] = u[x, Ny-1]
    
    # Scheme for all interior points
    u[1:-1,1:-1] = c1*(2*u_1[1:-1,1:-1] - c2*u_2[1:-1,1:-1] + \
      0.5*Cx2*((q[1:-1,1:-1] + q[2:,1:-1])*(u_1[2:,1:-1]-u_1[1:-1,1:-1]) - \
      (q[1:-1,1:-1] + q[:-2,1:-1])*(u_1[1:-1,1:-1] - u_1[:-2,1:-1])) + \
      0.5*Cy2*((q[1:-1,1:-1] + q[1:-1,2:])*(u_1[1:-1,2:] - u_1[1:-1,1:-1]) - \
      (q[1:-1,1:-1] + q[1:-1,:-2])*(u_1[1:-1,1:-1] - u_1[1:-1,:-2])) + dt2*f[1:-1,1:-1])
      
    if exact is not None:
        Err = exact[1:-1,1:-1] - u[1:-1,1:-1]  # Compare numerical soluion with exact solution
        E = sqrt(dx*dx*sum(Err**2))
    else:
        E = 0
    
    u_2 = u_1.copy()
    u_1 = u.copy()
    return u_1, u_2, E
