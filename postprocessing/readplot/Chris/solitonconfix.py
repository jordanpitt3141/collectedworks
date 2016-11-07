import csv
from numpy.linalg import norm
from scipy import *
from pylab import plot, show, legend,xlim,ylim,savefig,title,xlabel,ylabel,clf, loglog

order = "3"
num = "6"
wdir = "../../../../data/Joefemm/femacc/o3fem/sol/savenorms.txt"
sdir = "../../../results/Chrisfem/solcon/"



h0 = 0.1
g = 9.81
         
s = wdir
with open(s,'r') as file1:
    readfile = csv.reader(file1, delimiter = ',', quotechar='|', quoting=csv.QUOTE_MINIMAL)
        
    hL1 = []
    dxs = []
    uL1 = []
    ht = []
    ut = []
    j = -1
    for row in readfile:       
        if (j >= 0):
            dxs.append(float(row[0]))
            hL1.append(float(row[1]))
            uL1.append(float(row[2]))
                
        j = j + 1
    uL1 = array(uL1)
    hL1 = array(hL1)
    dxs = array(dxs)

a1 = 1.0
a0 = 10.0
g = 9.81
t0 = 0.0
n = len(dxs)  

s = sdir + "hL1.dat"
with open(s,'a') as file1:
    for i in range(0,n):
        s ="%3.8f%5s%1.15f\n" %(dxs[i]," ",hL1[i])
        file1.write(s)
s = sdir + "uL1.dat"
with open(s,'a') as file2:
    for i in range(0,n):
        s ="%3.8f%5s%1.15f\n" %(dxs[i]," ",uL1[i])
        file2.write(s)    