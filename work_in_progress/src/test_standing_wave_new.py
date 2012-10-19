from wave_solver import solver
from numpy import *
import math as m
import os, glob
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
    
C = 1
#Nx = 20
#Ny = 20
Lx = 4
Ly = 4
T = 4
c = 1.1
b = 0.05
w = pi
my = 2
mx = 2
version = "vectorized"
q_const = 0.4
cx = mx*pi/Lx
cy = my*pi/Ly
#dt = (1/c)*(1/(sqrt(1/0.04**2 + 1/0.04**2)))
dt = 0.0155


def exact_standing_wave(x,y,b,t,w,version="vectorized"):
    if version=="scalar":
        return m.exp(-b*t)*m.cos(mx*x*pi/Lx)*m.cos(my*y*pi/Ly)*m.cos(w*t)
    else:
        return m.exp(-b*t)*cos(mx*x*pi/Lx)*cos(my*y*pi/Ly)*m.cos(w*t)

        
    
def test_standing_wave(version, w, Nx):
    Ny = Nx
   
    def q(x,y):
        return ones((len(x),len(y))) * q_const
    
    def f(x,y,t):
        if version=="scalar":
            return m.exp(-b*t)*m.cos(cx*x)*m.cos(cy*y)*(m.cos(w*t)*(q_const*cx**2 + q_const*cy**2 - w**2) + b*w*m.sin(w*t))
        else:
           return exp(-b*t)*cos(cx*x)*cos(cy*y)*(cos(w*t)*(q_const*cx**2 + q_const*cy**2 - w**2) + b*w*sin(w*t))
        
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
        
    E_list, u, dx = solver(Lx, Ly, Nx, Ny, T, dt, c, I, q, V, f, b, version, w, exact_standing_wave, standing=True, make_plot=False)
    return E_list, dx
    

    
def compute_error(w, Nx):
    e_list, dx = test_standing_wave("vectorized", w, Nx)
    E = sqrt(dx*dx*sum(e_list**2))
    return E, dx
       

Nx_list = [20.0, 40.0, 80.0, 160.0]
Error_list=[]
E_dx_list=[]
dx_list=[]
for nx in Nx_list:
    error, dx = compute_error(w, nx)
    dx_list.append(dx)
    Error_list.append(error)
    E_dx_list.append(error/dx**2)
    #print "Nx = %g, E = %g" % (nx, error)
m = len(Nx_list)
print "dx-list:"
print dx_list
print "E/dx**2: "
print E_dx_list
print "Error_list:"
print Error_list
r = [log(Error_list[i-1]/Error_list[i])/log(dx_list[i-1]/dx_list[i]) for i in range(1, m, 1)]
noe = [Error_list[i-1]/Error_list[i] for i in range(1,m,1)]
print noe
print "Convergence rates: "
print r

    
    
def plot_exact():
    #Fy = 0.8
    #Fx = 0.8
    #Ft = 1/float(c) *1/sqrt(1/Fx**2 + 1/Fy**2)
    #dx = Fx*h
    #dy = Fy*h
    #dt = Ft*h
    #dt = (1/float(c))*(1/sqrt(1/dx**2 + 1/dy**2))
    N = int(round(float(T/dt)))
    x = linspace(0,Lx,Nx+1)
    y = linspace(0,Ly,Ny+1)
    dx = x[1] - x[0]
    dy = y[1] - y[0]
    X,Y = meshgrid(x,y)
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
    U_e0 = exact_standing_wave(X,Y,b,t[0],w)
    
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
        U_e = exact_standing_wave(X,Y,b,t[n],w)
        #print "t=%g" %n, U_e
        #print U_e
        savetxt("texttmp%.4d.txt"%n, U_e)

    
#plot_exact()