from wave_solver import solver
from numpy import *
import math as m
import os, glob
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
    
C = 1
Nx = 20
Ny = 20
Lx = 10
Ly = 10
T = 4
c = 1.1
b = 0.05
dt = -1.0
w = 1.3
my = 2
mx = 2
h = 0.1
version = "vectorized"
cx = mx*pi/Lx
cy = my*pi/Ly



def exact_manufactured_solution(x,y,b,t,w,version="vectorized"):
    if version=="scalar":
        return m.exp(-b*t)*m.cos(mx*x*pi/Lx)*m.cos(my*y*pi/Ly)*m.cos(w*t)
    else:
        return m.exp(-b*t)*cos(mx*x*pi/Lx)*cos(my*y*pi/Ly)*m.cos(w*t)
     
        
    
def test_manufactured_solution(version, h,w):
   
    def q(x,y):
            return x*y
    
    def f(x,y,t):
        if version=="scalar":
            f_val = m.exp(-b*t)*(m.cos(cx*x)*m.cos(cy*y)*(-w**2*m.cos(w*t) + b*w*m.sin(w*t)) + \
                cx*y*m.cos(cy*y)*m.cos(w*t)*(cx*x*m.cos(cx*x) + m.sin(cy*y)) + \
                cy*x*m.cos(cx*x)*m.cos(w*t)*(cy*y*m.cos(cy*y) + m.sin(cy*y)))
            return f_val
        else:
            f_val = exp(-b*t)*(cos(cx*x)*cos(cy*y)*(-w**2*cos(w*t) + b*w*sin(w*t)) + \
                cx*y*cos(cy*y)*cos(w*t)*(cx*x*cos(cx*x) + sin(cy*y)) + \
                cy*x*cos(cx*x)*cos(w*t)*(cy*y*cos(cy*y) + sin(cy*y)))
            return f_val
        
    def I(x,y):
        if version=="scalar":   
            return m.cos(cx*x)*m.cos(cy*y)
        else:
            return cos(cx*x)*cos(cy*y)
        
    def V(x,y):
        if version=="scalar":
            return -b*m.cos(cx*x)*m.cos(cy*y)
        else:
            return -b*cos(cx*x)*cos(cy*y)
        
    E_list, u = solver(Lx, Ly, Nx, Ny, T, dt, c, I, q, V, f, b, version, h, w, exact_manufactured_solution, standing=True, make_plot=False)
    return E_list
    

    
def compute_error(h, w):
    E_list = []
    e_list = test_manufactured_solution("vectorized", h, w)
    #print e_list[1]
    for e in e_list:
        temp = sqrt(h*sum(e**2))
        E_list.append(temp)
    s = 0
    for e in E_list:
        s += e**2
    E = sqrt(h*s)
    return E
  
  

#def test_omega():
   #wn = 10
   #error_test = 100
   #w_list = linspace(0.02,2*pi, wn+1)
   #for i in range(wn+1):
       #error = compute_error(h, w_list[i])
       #print "omega: %.4f" % w_list[i], "error: ", error
       #if error < error_test:
           #error_test = error
           #omega = w_list[i]
   #return error_test, omega
       
       

h_list = [8, 4, 2, 1, 0.5, 0.25]
Error_list=[]
for h_val in h_list:
    Error_list.append(compute_error(h_val, w))
print "E-list: ", Error_list
m = len(h_list)
r = [log(Error_list[i-1]/Error_list[i])/(log(h_list[i-1]/h_list[i])) for i in range(1, m, 1)]
print "Convergence rates: "
print r
    

    
    
"""
  m = len(dt_list)
    r = [log(E_values[i-1]/E_values[i])/
         log(dt_list[i-1]/dt_list[i])
         for i in range(1, m, 1)]
"""

    
    
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
            #U_e[i,j] = exact_manufactured_solution(x[i], y[j], b, t[0])
    U_e0 = exact_manufactured_solution(X,Y,b,t[0])
    
    #plt.savefig("tmp_%.4d.png" % 0)
    #ax.cla()
    savetxt("u0.txt", U_e0)
    for n in range(1,N+1):
        #    for i in range(Nx+1):
        #        for j in range(Ny+1):
        #            U_e[i,j] = exact_manufactured_solution(x[i], y[j], b, t[n])
        #ax.plot_surface(X,Y,U_e)
        #plt.savefig("tmp_%.4d.png" % n)
        #ax.cla()
        #print "Saved fig %g out of %g" % (n, N)
        U_e = exact_manufactured_solution(X,Y,b,t[n])
        #print "t=%g" %n, U_e
        #print U_e
        savetxt("texttmp%.4d.txt"%n, U_e)

        
#plot_exact()

