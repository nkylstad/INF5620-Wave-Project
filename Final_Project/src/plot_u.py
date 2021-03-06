import matplotlib.pyplot as plt
import scitools.easyviz as sci
import glob, os

def plot_u(u, x, t, n, b, Lx):
    plt.plot(x,u)
    plt.axis((0,Lx,-0.2,2.2))
    plt.title("t = %.4f,  b = %g" % (t, b))
    plt.savefig("tmp_%.4d.png" % n)
    print "Saving image %g" % n
    plt.clf()
    
def make_movie():
    sci.movie("tmp_*.png")
    for i in glob.glob("tmp_*.png"):
        os.remove(i)