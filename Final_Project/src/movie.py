from numpy import *
import math, os, subprocess, sys, glob
import matplotlib.pyplot as plt
from mayavi import mlab
import scitools.easyviz as sci

sci.movie("tmp_*.png")
for i in glob.glob("tmp_*.png"):
    os.remove(i)