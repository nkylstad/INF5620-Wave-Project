from mayavi import mlab
from numpy import *

a = zeros((51,51))
x = linspace(0,50,51)
y = linspace(0,50,51)
X,Y = meshgrid(x,y)

s = mlab.surf(a, X, Y, colormap = "bleh")