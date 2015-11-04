import csv
from numpy.linalg import norm
from scipy import *
from pylab import plot, show, legend,xlim,ylim,savefig,title,xlabel,ylabel,clf, loglog

order = "3"
num = "6"
wdir = "../../../results/Chrisfem/segurfem/scaledout20.txt"
sdir = "../../../results/Chrisfem/seg/"

gape = 2

h0 = 0.1
g = 9.81

        
s = wdir
with open(s,'r') as file1:
    readfile = csv.reader(file1, delimiter = ',', quotechar='|', quoting=csv.QUOTE_MINIMAL)
        
    h = []
    t = []
    j = -1
    for row in readfile:       
        if (j >= 0):
            dx =float(row[0])
            dt =float(row[1])
            t.append(float(row[2]))
            h.append(float(row[3]))
                
        j = j + 1
    t = array(t)
    h = 1.5*array(h)
#1.5 factor very important.
 
n = len(h)

s = sdir + "20.dat"
with open(s,'a') as file1:
    for i in range(0,n,gape):
        s ="%3.8f%5s%1.15f\n" %(t[i]," ",h[i])
        file1.write(s)    