import glob, os

#for i in glob.glob("tmp_*.png"):
    #os.remove(i)
    
os.remove("initial.txt")
os.remove("u0.txt")
for i in glob.glob("texttmp*.txt"):
    os.remove(i)
for j in glob.glob("wtmp*.png"):
    os.remove(j)