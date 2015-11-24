import csv
from numpy.linalg import norm
from scipy import *
from pylab import plot, show, legend,xlim,ylim,savefig,title,xlabel,ylabel,clf, loglog

order = "3"
num = "6"
wdir = "../../../../data/raw/dbchrisworking/o3femh11p8h01/outlast.txt"
sdir = "../../../../data/postprocessing/Chris/1p8and1/"

gape = 1

h0 = 0.1
g = 9.81
         
s = wdir
with open(s,'r') as file1:
    readfile = csv.reader(file1, delimiter = ',', quotechar='|', quoting=csv.QUOTE_MINIMAL)
        
    h = []
    x = []
    u = []
    ht = []
    ut = []
    j = -1
    for row in readfile:       
        if (j >= 0):
            dx =float(row[0])
            dt =float(row[1])
            t =float(row[2])
            x.append(float(row[3]))
            h.append(float(row[4]))
            u.append(float(row[6]))
                
        j = j + 1
    u = array(u)
    h = array(h)  
    
    
n = len(h)


if not os.path.exists(sdir):
    os.makedirs(sdir)
    
s = sdir + "h.dat"
with open(s,'a') as file1:
    for i in range(0,n,gape):
        s ="%3.8f%5s%1.15f\n" %(x[i]," ",h[i])
        file1.write(s)
s = sdir + "u.dat"
with open(s,'a') as file2:
    for i in range(0,n,gape):
        s ="%3.8f%5s%1.15f\n" %(x[i]," ",u[i])
        file2.write(s)  