import csv
from numpy.linalg import norm
from scipy import *
from pylab import plot, show, legend,xlim,ylim,savefig,title,xlabel,ylabel,clf, loglog
from matplotlib2tikz import save as tikz_save

#wdir = "../../../../data/raw/DBSWanew100s/o3/10/10.0/"
wdir = "../../../../data/raw/trackleadsola10new/o3/"
sdir = "../../../../data/postprocessing/DBSWfront/o3/10/10.0/"

if not os.path.exists(sdir):
    os.makedirs(sdir)

#wdir = "../../../data/Joe/dbh/o"+order+"/"
#sdir = "../../../../written/exportpic/dambreak/ex/"

#gap = 50
gap = 10

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
            x.append(float(row[4]))
            h.append(float(row[5]))
            u.append(float(row[7]))
                
        j = j + 1
    x = array(x)
    u = array(u)
    h = array(h)

xbeg = int((0)/dx)
xend = int((1000)/dx)
xt =x[xbeg:xend:gap]
ht =h[xbeg:xend:gap]


"""
n = len(xt)
s = sdir + "hDBSW.dat"
with open(s,'w') as file1:
    for i in range(n):
        s ="%3.8f%5s%1.15f\n" %(xt[i] + 200," ",ht[i])
        file1.write(s)
"""  


n = len(xt)
s = sdir + "hDB.dat"
with open(s,'w') as file1:
    for i in range(n):
        s ="%3.8f%5s%1.15f\n" %(xt[i]," ",ht[i])
        file1.write(s)
 