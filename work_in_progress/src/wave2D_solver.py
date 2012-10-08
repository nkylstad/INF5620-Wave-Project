from numpy import linspace, zeros
from math import *
from scitools.easyviz import *
import sys


def solver(I, V, q, Lx, Ly, Nx, Ny, c, dt):
    
    x = linspace(0,Lx,Nx+1)  # mesh points in x-direction
    y = linspace(0,Ly,Ny+1)  # mesh points in y-direction
    dx = x[0] - x[1]
    dy = y[0] - y[1]
    
    stability_limit = (1/float(c))*(1/sqrt(1/dx**2 + 1/dy**2))
    if dt <= 0:
        dt = 1*stability_limit
    elif dt > stability_limit:
        print 'error: dt=%g exceeds the stability limit %g' % \
              (dt, stability_limit)


    N = int(round((T/float(dt))))
    t = linspace(0,N*dt, N+1)  # mesh points in time
    Cx2 = (dt/(sqrt(2)*dx))**2; Cy2 = (dt/(sqrt(2)*dy))**2
    c1 = 1/(1 + ((b*dt)/2))
    c2 = 1 + (b*dt)/2
    
    up = zeros(Nx+1,Ny+1)
    u = zeros(Nx+1, Ny+1)
    um = zeros(Nx+1, Ny+1)
    
    # Initial conditions
    for i in range(Nx+1):
        for j in range(Ny+1):
            u[i][j] = I(x[i],y[j])
            
    # First step:
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
            
            up[i][j] = c1*(2*u[i][j] - c2*(u[i][j] - dt*V(x[i],y[j])) + Cx2*((q(x[i],y[j]) + q(x[spx],y[j]))*\
            (u[spx][j] - u[i][j]) - (q(x[smx],y[j]) + q(x[i],y[j]))*(u[i][j] - u[smx][j])) + \
            Cy2*((q(x[i],y[j]) + q(x[i],y[spy]))*(u[i][spy] - u[i][j]) - (q(x[i],y[j]) + q(x[i],y[smy]))*\
            (u[i][j] - u[i][smy])) + f(x[i],y[j],t[0]))

print "No syntax errors"

   