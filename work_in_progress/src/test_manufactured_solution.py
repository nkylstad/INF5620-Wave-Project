from wave_solver import solver
from numpy import *
import math as m
import os, glob


#os.remove("initial.txt")
#os.remove("u0.txt")
#for i in glob.glob("texttmp*.txt"):
    #os.remove(i)
#for j in glob.glob("wtmp*.png"):
    #os.remove(j)

    
C = 1
Nx = 50
Ny = 50
Lx = 4
Ly = 4
T = 10
c = 1.1
b = 0.3
dt = -1.0
w = 1
my = 2
mx = 2
#h = 0.1
version = "vectorized"


def exact_manufactured_solution(x,y,b,t,version="vectorized"):
    if version=="scalar":
        return m.exp(-b*t)*m.cos(mx*x*pi/Lx)*m.cos(my*y*pi/Ly)*m.cos(w*t)
    else:
        return exp(-b*t)*cos(mx*x*pi/Lx)*cos(my*y*pi/Ly)*cos(w*t)
        #return 
        
    
def test_manufactured_solution(version, h):
   
    q_const = ((w*Lx*Ly)/pi)**2*(1/((Lx*my)**2 + (Ly*mx)**2))
    q = ones((Nx+1,Ny+1)) * q_const
    
    def f(x,y,t):
        if version=="scalar":
            return -w*b*m.exp(-b*t)*m.cos(mx*x*pi/Lx)*m.cos(my*y*pi/Ly)*m.sin(w*t)
        else:
            
            return -w*b*exp(-b*t)*cos(mx*x*pi/Lx)*cos(my*y*pi/Ly)*sin(w*t)
        
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
        
    
    E_list = solver(Lx, Ly, Nx, Ny, T, dt, c, I, q, V, f, b, version, h, exact_manufactured_solution, standing=True, make_plot=False)
    return E_list
    
def compute_error(h):
    E_list = []
    e_list = test_manufactured_solution("vectorized", h)
    for e in e_list:
        temp = sqrt(sum(e**2))
        E_list.append(temp)
    s = 0
    print E_list
    for e in E_list:
        s += h*e**2
    E = sqrt(s)
    return E/(h**2)
    
#h_list = [1.0, 0.5, 0.1, 0.05, 0.01]
h_list = [0.5]
#for h_val in h_list:
    #error = compute_error(h_val)
    #print "h = %g, E/h^2 = %g" % (h_val, error)

        
    

def plot_exact():
    x = linspace(0,Lx,Nx+1)
    y = linspace(0,Ly,Ny+1)
    X,Y = meshgrid(x,y)
    dx = x[1]-x[0]
    dy = y[1]-y[0]
    dt = 0.05 #(1/float(c))*(1/sqrt(1/dx**2 + 1/dy**2))
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
    U_e0 = exact_manufactured_solution(X,Y,b,t[0])
    savetxt("u0.txt", U_e0)
    for n in range(1,N+1):
        U_e = exact_manufactured_solution(X,Y,b,t[n])
        #print U_e
        savetxt("texttmp%.4d.txt"%n, U_e)

        
plot_exact()
