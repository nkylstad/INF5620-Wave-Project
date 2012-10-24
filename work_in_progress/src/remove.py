import glob, os

# for i in glob.glob("tmp_*.png"):
#     os.remove(i)
    
os.remove("initial.txt")
os.remove("u0.txt")
os.remove("hill.txt")
os.remove("H0.txt")
os.remove("q.txt")
for i in glob.glob("texttmp*.txt"):
   os.remove(i)
for j in glob.glob("u*.png"):
   os.remove(j)