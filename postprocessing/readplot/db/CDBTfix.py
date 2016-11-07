import csv
import os
from numpy.linalg import norm
from scipy import *
from pylab import plot, show, legend,xlim,ylim,savefig,title,xlabel,ylabel,clf, loglog

wdirbase = "../../../../data/raw/CDBThetaHAMILYTIME/o2/"
sdirbase = "../../../../data/postprocessing/CDBCNF/"

Evals = []
thetas = []
for i in range(0,11):
#for i in range(1,2):
    #k = 1 + i*0.5
    k = 1 + i* 0.1
    wdir = wdirbase + str(k) + "/"    
    sdir = sdirbase + "theta="+str(k) + "/"   
    
    if not os.path.exists(sdir):
        os.makedirs(sdir)
    
             
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
                Evali =float(row[3])
                Evalf =float(row[4])
                x.append(float(row[5]))
                h.append(float(row[6]))
                u.append(float(row[8]))
                    
            j = j + 1
        x = array(x)
        u = array(u)
        h = array(h)
    
    relEval = abs(Evali - Evalf)/ abs(Evali)
    Evals.append(relEval)
    thetas.append(k)
    
    n = len(x)
    s = sdir + "XAndH.dat"
    with open(s,'w') as file1:
        for i in range(n):
            s ="%3.8f%5s%1.15f\n" %(x[i]," ",h[i])
            file1.write(s)
            
    s = sdir + "XAndU.dat"
    with open(s,'w') as file1:
        for i in range(n):
            s ="%3.8f%5s%1.15f\n" %(x[i]," ",u[i])
            file1.write(s)
s = sdirbase + "ThetaAndHamiltonianConservation.dat"
n = len(thetas)
with open(s,'w') as file1:
    for i in range(n):
        s ="%3.8f%5s%1.15f\n" %(thetas[i]," ",Evals[i])
        file1.write(s)  