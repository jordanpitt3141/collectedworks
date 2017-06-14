import csv
from numpy.linalg import norm
from scipy import *
from pylab import plot, show, legend,xlim,ylim,savefig,title,xlabel,ylabel,clf, loglog
import os
wdir = "../../../../data/raw/Joebigsmooth/o3/11/12/"

wdir2 = "../../../../data/raw/Joebigsmooth/o3/11/20/"
sdir = "../../../../data/postprocessing/PRES/DBstruct/12/"

if not os.path.exists(sdir):
    os.makedirs(sdir)

#wdir = "../../../data/Joe/dbh/o"+order+"/"
#sdir = "../../../../written/exportpic/dambreak/ex/"

#gap = 50
gap = 8

h0 = 0.1
g = 9.81

         
s = wdir + "outlast.txt"
with open(s,'r') as file1:
    readfile = csv.reader(file1, delimiter = ',', quotechar='|', quoting=csv.QUOTE_MINIMAL)
        
    h = []
    u = []
    x = []
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
    x = array(x)
    u = array(u)
    h = array(h)

xbeg = int((300)/dx)
xend = int((700)/dx)
xt =x[xbeg:xend:gap]
ht =h[xbeg:xend:gap]
"""
s = wdir2 + "outlast.txt"
with open(s,'r') as file1:
    readfile = csv.reader(file1, delimiter = ',', quotechar='|', quoting=csv.QUOTE_MINIMAL)
        
    h = []
    u = []
    x = []
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
    x = array(x)
    u = array(u)
    h = array(h)

xbeg = int((350)/dx)
xend = int((650)/dx)
xt1 =x[xbeg:xend:gap]
ht1 =h[xbeg:xend:gap]
"""
"""
plot(xt,ht ,'-b')
xlim([350,650])
ylim([1.0,1.8])
xlabel("$x$ ($m$)")
ylabel("$h$ ($m$)")
"""
#s = sdir + "o" + order +"n" + num +".tikz" 
#tikz_save(s);      
#clf();

n = len(xt)
s = sdir + "hz.dat"
with open(s,'w') as file1:
    for i in range(n):
        s ="%3.8f%5s%1.15f\n" %(xt[i]," ",ht[i])
        file1.write(s)
  