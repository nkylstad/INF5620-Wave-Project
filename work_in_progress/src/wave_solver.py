from numpy import *
import math, os, subprocess, sys, glob
import matplotlib.pyplot as plt
from mayavi import mlab
import scitools.easyviz as sci

def solver(Lx, Ly, Nx, Ny, T, dt, c, I, q, V, f, b, version, oneD=False):
    
    x = linspace(0,Lx,Nx+1)
    y = linspace(0,Ly,Ny+1)
    X,Y = meshgrid(x,y)
    dx = float(x[1] - x[0]) 
    dy = float(y[1] - y[0])
   
    stability_limit = (1/float(c))*(1/sqrt(1/dx**2 + 1/dy**2))
    if dt <= 0:
        dt = 1*stability_limit
    elif dt < stability_limit:
        print "Error: dt too large."
    
    if oneD:
        y = zeros(Ny+1)
        dy = 1
        Cy2 = 0
        dt = dx/c
        Y = meshgrid(y, y)
        
    N = int(round(float(T/dt)))
    t = linspace(0,T,N+1)
    c1 = 1/(1 + (b*dt)/2)
    c2 = 1 - (b*dt)/2
    Cx2 = (dt/dx)**2
    Cy2 = (dt/dy)**2
    dt2 = dt**2
    
    u = zeros((Nx+1, Ny+1))  # The new soluion at the next timestep
    u_1 = zeros((Nx+1, Ny+1)) # The solution from the current time step
    u_2 = zeros((Nx+1, Ny+1))  # The solution from the previous time step
    
    #s = mlab.mesh(X[1:-1,1:-1], Y[1:-1,1:-1], u[1:-1,1:-1])
    # Initial conditions
    if version == "scalar":
        for i in range(0,Nx+1):
            for j in range(0,Ny+1):
                u_2[i,j] = I(x[i],y[j])
    else:  # vectorized version
        u_2[:,:] = I(X,Y)
        
    plt.plot(x,u_2)
    plt.axis((x[0],x[-1],-0.2,2.2))
    plt.savefig("tmp_%.4d.png" % 0)
    plt.clf()
        
        
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
        u_1[1:-1,1:-1] = u_2[1:-1,1:-1] + dt*c2*V[1:-1,1:-1] + \
            0.5**2*Cx2*((q[2:,1:-1] + q[1:-1,1:-1])*(u_2[2:,1:-1] - u_2[1:-1,1:-1]) - (q[1:-1,1:-1] + q[:-2,1:-1])*(u_2[1:-1,1:-1] - u_2[:-2,1:-1])) + \
            0.5**2*Cy2*((q[1:-1,2:] + q[1:-1,1:-1])*(u_2[1:-1,2:] - u_2[1:-1,1:-1]) - (q[1:-1,1:-1] + q[1:-1,:-2])*(u_2[1:-1,1:-1] - u_2[1:-1,:-2]))
            # + 0.5*dt2*f(xkjbsdjkzcv)
            
        # boundary conditions:
        u_1[0,:] = u_1[1,:]
        u_1[-1,:] = u_1[-2,:]
        u_1[:,0] = u_1[:,1]
        u_1[:,-1] = u_1[:,-2]
        
    plt.plot(x,u_1)
    plt.axis((x[0],x[-1],-0.2,2.2))
    plt.savefig("tmp_%.4d.png" % 1)
    plt.clf()
    E = zeros(N+1)

    for n in range(2,N+1):
        
        if version == "scalar":
            u_1,u_2, E_val = advance_scalar(u, u_1, u_2, Nx, Ny, x, y, q, f, c1, c2, Cx2, Cy2, dt2, b, t[n], exact_manufactured_solution)
            #plt.figure()
            plt.plot(x,u_1)
            plt.show()
            plt.clf()
                 
        else:
            u_1,u_2 = advance_vectorized(u, u_1, u_2, q, f, c1, c2, Cx2, Cy2, t, exact_manufactured_solution)
            #plt.figure()
            
            plt.plot(x,u_1)
            plt.axis((x[0],x[-1],-0.2,2.2))
            plt.savefig("tmp_%.4d.png" % n)
            plt.show()
            plt.clf()
            #plot_u(u_1, x, X, y, Y, t, n, 2)
    sci.movie("tmp_*.png")
    for i in glob.glob("tmp_*.png"):
        os.remove(i)
    return u_1       
         
         
def advance_scalar(u, u_1, u_2, Nx, Ny, x, y, q, f, c1, c2, Cx2, Cy2, dt2, b, tn, exact):
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
                
            u_exact = exact(x[i],y[j],b,tn)
            temp = abs(u_exact - u[i,j])
            if temp > error:
                error = temp
            
    E_val = error
                
    u_2 = u_1.copy()
    u_1 = u.copy()
    return u_1, u_2 , E_val
    
def advance_vectorized(u, u_1, u_2, q, f, c1, c2, Cx2, Cy2, t, exact_manufactured_solution):
    u[1:-1,1:-1] = c1*(2*u_1[1:-1,1:-1] - c2*u_2[1:-1,1:-1] + \
      0.5*Cx2*((q[1:-1,1:-1] + q[2:,1:-1])*(u_1[2:,1:-1]-u_1[1:-1,1:-1]) - \
      (q[1:-1,1:-1] + q[:-2,1:-1])*(u_1[1:-1,1:-1] - u_1[:-2,1:-1])) + \
      0.5*Cy2*((q[1:-1,1:-1] + q[1:-1,2:])*(u_1[1:-1,2:] - u_1[1:-1,1:-1]) - \
      (q[1:-1,1:-1] + q[1:-1,:-2])*(u_1[1:-1,1:-1] - u_1[1:-1,:-2])) + f[1:-1,1:-1])
      
    #boundary conditions:
    
    u[0,:] = u[1,:]
    u[-1,:] = u[-2,:]
    u[:,0] = u[:,1]
    u[:,-1] = u[:,-2]
    
    u_2 = u_1.copy()
    u_1 = u.copy()
    return u_1, u_2
    
    
def plot_u(u, x, X, y, Y, t, n, plot_method):
    """if t[n] == 0:
    time.sleep(2)
    if plot_method == 1:
    mesh(x, y, u, title='t=%g' % t[n], zlim=[-1,1],
    caxis=[-1,1])
    elif plot_method == 2:
    surfc(xv, yv, u, title='t=%g' % t[n], zlim=[-1, 1],
    colorbar=True, colormap=hot(), caxis=[-1,1],
    shading='flat')
    
    fig = plt.figure()
    ax = fig.add_subplot(111, projection="3d")
    Axes3D.plot_surface(X,Y,u)
    if plot_method > 0:
        time.sleep(0) # pause between frames
        filename = 'tmp_%04d.png' % n
        #savefig(filename)  # time consuming - dropped
    """
"""
    def plot_3D(x, y, u, s):
    #mlab.mesh(x,y,u)
    #mlab.surf(x,y,u)
    s.mlab_source.scalars = u[1:-1,1:-1]
    #mlab.savefig("wtmp%.5f.png" % n)
"""
    
pi = math.pi
w = 0.8
Lx = 5
Ly = 5
Nx = 20
Ny = 20
T = 4
my = 2
mx = 2
b = 0.1
c = 1.1
dt = -1    
    
def exact_manufactured_solution(x,y,b,t):
    return exp(-b*t)*cos(mx*x*pi/Lx)*cos(my*y*pi/Ly)*sin(w*t)
    
def test_manufactured_solution(version):

    
    q_const = ((w*Lx*Ly)/pi)**2*(1/((Lx*my)**2 + (Ly*mx)**2))
    q = ones((Nx+1,Ny+1)) * q_const
    
    def f(x,y,t):
        return w*b*exp(-b*t)*cos(mx*x*pi/Lx)*cos(my*y*pi/Ly)*sin(w*t)
        
    def I(x,y):
        return cos(mx*x*pi/Lx)*cos(my*y*pi/Ly)
        
    def V(x,y):
        return -b*cos(mx*x*pi/Lx)*cos(my*y*pi/Ly)
        
    
    E = solver(Lx, Ly, Nx, Ny, T, dt, c, I, q, V, f, b, version)
    
    print E

    
#test_manufactured_solution("scalar")


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
        
    
    u = solver(Lx,Ly,Nx,Ny,T,dt,c,I,q,V,f,b,version)
    print "Final solution!!!", u
 
#test_constant("scalar")

def plot_1D(u, x, fig, ax, n):
    """
    ax.cla()
    ax.imshow(rand(5,5)interpolation="nearest")
    fname = "tmp1D_%.4f.png" % n
    return fname
    """
    plt.plot(x,u,'b')
    plt.axis((x[0],x[-1],-0.5,3))
    plt.xlabel('x')
    plt.ylabel('u')
    plt.title('1D plug wave', fontsize=20)

    fname = "tmp1D_%.4f.png" % n 
    return fname

    
def test_1D_plug(version):
    C = 1
    Nx = 200
    Ny = 200
    Lx = 4
    Ly = 4
    T = 20
    c = 1.
    b =0
    dt = 0.04
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
test_1D_plug("vectorized")

    
