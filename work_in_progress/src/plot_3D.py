from numpy import *
import math as m
from mayavi import mlab
from scitools.easyviz import *

f = open("initial.txt", 'r')
Nx=0; Ny=0; Lx=0; Ly=0; T=0; N=0
for line in f:
    words = line.split("=")
    if words[0]=="Nx":
        Nx = int(words[1])
    elif words[0]=="Ny":
        Ny = int(words[1])
    elif words[0]=="Lx":
        Lx = int(words[1])
    elif words[0]=="Ly":
        Ly = int(words[1])
    elif words[0]=="T":
        T = int(words[1])
    elif words[0]=="N":
        N = int(words[1])
    else:
        print "Error in file initial.txt"
print "Finished with initial conditions"
f.close()

f = open("H0.txt", 'r')
H0 = float(f.readline())
f.close()

x = linspace(0,Lx,Nx+1)
y = linspace(0,Ly,Ny+1)
X,Y = meshgrid(x,y)
t = linspace(0,T,N+1)

f = open("hill.txt", 'r')
hill = zeros((Nx+1,Ny+1))
i = 0
for line in f:
    numbers = line.split()
    hill[i,:] = array([float(numbers[j])-H0 for j in range(len(numbers))])
    i += 1
f.close()

f = open("u0.txt", 'r')
u0 = zeros((Nx+1, Ny+1))
i = 0
for line in f:
    numbers = line.split()
    u0[i,:] = array([float(numbers[j]) for j in range(len(numbers))])
    i += 1
f.close()

print "Plotting image 1 of %g" % N
s = mlab.mesh(X, Y, u0,colormap='winter', opacity=0.8)  # color=(0.1,0.7,0.8)
hill_p = mlab.mesh(X,Y,hill, colormap='gray', opacity=0.8)
mlab.title("t=%.4f" % t[0], color=(1,1,1),size=0.5)
mlab.options.offscreen = True
mlab.savefig("u%.3d.png" % 0)
mlab.clf()

for n in range(1,N):
    if n%4==0:
        filename = "texttmp%.4d.txt" % n
        f = open(filename, 'r')
        u = zeros((Nx+1, Ny+1))
        i = 0
        for line in f:
            numbers = line.split()
            u[i,:] = array([float(numbers[j]) for j in range(len(numbers))])
            i += 1
        f.close()
        print "Plotting image %g of %g" % (n+1, N)
        s = mlab.mesh(X,Y,u,colormap='winter', opacity=0.8)
        hill_p = mlab.mesh(X,Y,hill, colormap='gray', opacity=0.7)
        mlab.title("t=%.4f" % t[n], color=(1,1,1),size=0.5)
        mlab.options.offscreen = True
        mlab.savefig("u%.4d.png" % n)
        mlab.clf()
    
import os, glob

movie("u*.png", fps=10)
for i in glob.glob("u*.png"):
    os.remove(i)
