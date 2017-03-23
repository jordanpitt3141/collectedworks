import csv
from numpy.linalg import norm
from scipy import *
from pylab import plot, show, legend,xlim,ylim,savefig,title,xlabel,ylabel,clf, loglog
import os

time = 10

wdir = "../../../../data/raw/solbumpex/o2/16/"
sdir = "../../../../data/postprocessing/solbumpex/o2/40s/"


if not os.path.exists(sdir):
    os.makedirs(sdir)

s = wdir + "out345176.txt"
with open(s,'r') as file1:
    readfile = csv.reader(file1, delimiter = ',', quotechar='|', quoting=csv.QUOTE_MINIMAL)
        
    h = []
    u = []
    b = []
    x = []
    j = -1
    for row in readfile:       
        if (j >= 0):
            dx =float(row[0])
            dt =float(row[1])
            t =float(row[2])
            Evaln = float(row[3])
            x.append(float(row[4]))
            h.append(float(row[5]))
            u.append(float(row[7]))
            b.append(float(row[8]))
               
        j = j + 1

h = array(h)
b = array(b)
x = array(x)
u = array(u)


n = len(x)
s = sdir + "stage.dat"
with open(s,'w') as file1:
    for i in range(n):
        s ="%3.8f%5s%1.15f\n" %(x[i]," ",h[i] + b[i])
        file1.write(s)
        
s = sdir + "h.dat"
with open(s,'w') as file1:
    for i in range(n):
        s ="%3.8f%5s%1.15f\n" %(x[i]," ",h[i])
        file1.write(s)
s = sdir + "u.dat"
with open(s,'w') as file2:
    for i in range(n):
        s ="%3.8f%5s%1.15f\n" %(x[i]," ",u[i])
        file2.write(s)

s = sdir + "bed.dat"
with open(s,'w') as file1:
    for i in range(n):
        s ="%3.8f%5s%1.15f\n" %(x[i]," ",b[i])
        file1.write(s)

