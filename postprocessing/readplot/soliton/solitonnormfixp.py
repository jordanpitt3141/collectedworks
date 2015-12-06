import csv
from numpy.linalg import norm
from scipy import *
from pylab import plot, show, legend,xlim,ylim,savefig,title,xlabel,ylabel,clf, loglog

#wdir = "../../../../data/JoeFDaccfin/FDacc/"
#sdir = "../../../results/FDsolacc/FD/"

wdir = "../../../../data/raw/nsolcon/o2af/"
sdir = "../../../../data/postprocessing/nsolcon30/o2af/"

if not os.path.exists(sdir):
    os.makedirs(sdir)

h0 = 0.1
g = 9.81
         
s = wdir + "savenorms.txt"
with open(s,'r') as file1:
    readfile = csv.reader(file1, delimiter = ',', quotechar='|', quoting=csv.QUOTE_MINIMAL)
        
    h = []
    u = []
    dxs = []
    j = -1
    for row in readfile:       
        if (j >= 0):
            dxs.append(float(row[0]))
            h.append(float(row[1]))
            u.append(float(row[2]))
                
        j = j + 1
    u = array(u)
    h = array(h)
    dxs = array(dxs)
    
         
#scaling
ldx = dxs
lh = h
lu = u

    
n = len(lh)

s = sdir + "h.dat"
with open(s,'a') as file1:
    for i in range(n):
        s ="%3.8f%5s%1.15f\n" %(ldx[i]," ",lh[i])
        file1.write(s)
s = sdir + "u.dat"
with open(s,'a') as file2:
    for i in range(n):
        s ="%3.8f%5s%1.15f\n" %(ldx[i]," ",lu[i])
        file2.write(s)
