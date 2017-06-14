import csv
from numpy.linalg import norm
from scipy import *
from pylab import plot, show, legend,xlim,ylim,savefig,title,xlabel,ylabel,clf, loglog
from numpy import zeros
import os

wdir = "../../../../data/raw/trackleadsol/o3/"
sdir = "../../../../data/postprocessing/RLAplusRel/o3/"

if not os.path.exists(sdir):
    os.makedirs(sdir)

#wdir = "../../../data/Joe/dbh/o"+order+"/"
#sdir = "../../../../written/exportpic/dambreak/ex/"

#gap = 50
gap = 1

h0 = 0.1

         
s = wdir + "aplus.txt"
with open(s,'r') as file1:
    readfile = csv.reader(file1, delimiter = ',', quotechar='|', quoting=csv.QUOTE_MINIMAL)
        
    a = []
    x = []
    t = []
    g = []
    j = -1
    for row in readfile:       
        if (j >= 0):
            x.append(float(row[0]))
            t.append(float(row[1]))
            a.append(float(row[2]))
            g.append(1.7399758)
                
        j = j + 1
    x = array(x)
    t = array(t)
    a = array(a)
    g = array(g)

n = len(x)
c = zeros(n)    
for i in range(1,n):
    c[i] =  (x[i] - x[i-1]) / (t[i] - t[i-1])

s = sdir + "a.dat"
with open(s,'w') as file1:
    for i in range(n):
        s ="%3.8f%5s%1.15f\n" %(t[i]," ",(a[i] - g[i])/g[i] )
        file1.write(s)
        
s = sdir + "s.dat"
with open(s,'w') as file1:
    for i in range(n):
        s ="%3.8f%5s%1.15f\n" %(t[i]," ",c[i])
        file1.write(s)
        
s = sdir + "grim.dat"
with open(s,'w') as file1:
        s ="%3.8f\n" %(g[i])
        file1.write(s)
  