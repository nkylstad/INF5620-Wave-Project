from wave_solver import solver
from numpy import *
import math as m
import os, glob
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
    
C = 1
Nx = 20
Ny = 20
Lx = 4
Ly = 4
T = 10
c = 1.1
b = 0.05
dt = -1.0
w = 1.1
my = 2
mx = 2
h = 0.1
version = "vectorized"


def exact_standing_wave(x,y,b,t,w,version="vectorized"):
    if version=="scalar":
        return m.exp(-b*t)*m.cos(mx*x*pi/Lx)*m.cos(my*y*pi/Ly)*m.cos(w*t)
    else:
        return m.exp(-b*t)*cos(mx*x*pi/Lx)*cos(my*y*pi/Ly)*m.cos(w*t)

        
    
def test_standing_wave(version, h, w):
   
    def q(x,y):
        q_const = ((w*Lx*Ly)/pi)**2*(1/((Lx*my)**2 + (Ly*mx)**2))
        return ones((len(x),len(y))) * q_const
    
    def f(x,y,t):
        if version=="scalar":
            return w*b*m.exp(-b*t)*m.cos(mx*x*pi/Lx)*m.cos(my*y*pi/Ly)*m.sin(w*t)
        else:
            if t==0:
                return zeros((len(x),len(y)))
            else:
                return w*b*exp(-b*t)*cos(mx*x*pi/Lx)*cos(my*y*pi/Ly)*sin(w*t)
        
    def I(x,y):
        if version=="scalar":   
            return m.cos(mx*x*pi/Lx)*m.cos(my*y*pi/Ly)
        else:
            return cos(mx*x*pi/Lx)*cos(my*y*pi/Ly)
        
    def V(x,y):
        if version=="scalar":
            return -b*m.cos(mx*x*pi/Lx)*m.cos(my*y*pi/Ly)
        else:
            return -b*cos(mx*x*pi/Lx)*cos(my*y*pi/Ly)
        
    E_list, u = solver(Lx, Ly, Nx, Ny, T, dt, c, I, q, V, f, b, version, h, w, exact_standing_wave, standing=True, make_plot=False)
    return E_list
    

    
def compute_error(h, w):
    e_list = test_standing_wave("vectorized", h, w)
    E = sqrt(h*sum(e_list**2))
    return E/(h**2)
    #return E
       
       

h_list = [2, 1, 0.5, 0.25, 0.125]
for h_val in h_list:
    error = compute_error(h_val, w)
    print "h = %g, E/h^2 = %g" % (h_val, error)
#Error_list=[]
#for h_val in h_list:
    #Error_list.append(compute_error(h_val, w))
#print "E-list: ", Error_list
#m = len(h_list)
#r = [log(Error_list[i-1]/Error_list[i])/(log(h_list[i-1]/h_list[i])) for i in range(1, m, 1)]
#print "Convergence rates: "
#print r

    
    
def plot_exact():
    x = linspace(0,Lx,Nx+1)
    y = linspace(0,Ly,Ny+1)
    X,Y = meshgrid(x,y)
    Fx = 0.8
    Fy = 0.8
    Ft = 1/c *1/sqrt(1/Fx**2 + 1/Fy**2)
    dx = Fx*h
    dy = Fy*h
    dt = Ft*h
    #dt = (1/float(c))*(1/sqrt(1/dx**2 + 1/dy**2))
    N = int(round(float(T/dt)))
    t = linspace(0,T,N+1)
    file = open("initial.txt",'w')
    file.write("Nx="+str(Nx)+"\n")
    file.write("Lx="+str(Lx)+"\n")
    file.write("Ny="+str(Ny)+"\n")
    file.write("Ly="+str(Ly)+"\n")
    file.write("T="+str(T)+"\n")
    file.write("N="+str(N)+"\n")
    file.close()
    U_e = zeros((Nx+1,Ny+1))
    #fig = plt.figure(111)
    #fig.clf()
    #ax = fig.gca(projection='3d')
    #for i in range(Nx+1):
        #for j in range(Ny+1):
            #U_e[i,j] = exact_standing_wave(x[i], y[j], b, t[0])
    U_e0 = exact_standing_wave(X,Y,b,t[0])
    
    #plt.savefig("tmp_%.4d.png" % 0)
    #ax.cla()
    savetxt("u0.txt", U_e0)
    for n in range(1,N+1):
        #    for i in range(Nx+1):
        #        for j in range(Ny+1):
        #            U_e[i,j] = exact_standing_wave(x[i], y[j], b, t[n])
        #ax.plot_surface(X,Y,U_e)
        #plt.savefig("tmp_%.4d.png" % n)
        #ax.cla()
        #print "Saved fig %g out of %g" % (n, N)
        U_e = exact_standing_wave(X,Y,b,t[n])
        #print "t=%g" %n, U_e
        #print U_e
        savetxt("texttmp%.4d.txt"%n, U_e)

        
#plot_exact()

