import csv
from numpy.linalg import norm
from scipy import *
from pylab import plot, show, legend,xlim,ylim,savefig,title,xlabel,ylabel,clf, loglog
import os

wdirbase = "../../../../data/raw/CDBThetaHAMILYTIME/o2/"
sdirbase = "../../../../data/postprocessing/CDBThetaHamilOverTime/dx=0.012/"


thetas = [0,1.0,1.1,1.2,1.3,1.4,1.5,1.6,1.7,1.8,1.9,2.0]
#thetas = [2.0]
Evals = []
#for i in range(0,11):
for theta in thetas:
    wdir = wdirbase + str(theta) + "/"    
    sdir = sdirbase + "theta="+str(float(theta)) + "/"   
    
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
        
                 
    s = wdir + "Hamil.txt"
    with open(s,'r') as file2:
        readfile = csv.reader(file2, delimiter = ',', quotechar='|', quoting=csv.QUOTE_MINIMAL)
            
        times = []
        Hamils = []
        j = -1
        for row in readfile:       
            if (j >= 0):
                dx =float(row[0])
                thetaF =float(row[1])
                times.append(float(row[2]))
                Hamils.append(float(row[3]))
                    
            j = j + 1
        times = array(times)
        Hamils = array(Hamils)
    
    relEval = abs(Evali - Evalf)/ abs(Evali)
    Evals.append(relEval)
    
    m = len(times)        
    s = sdir + "TAndHamiltonian.dat"
    with open(s,'w') as file1:
        for i in range(m):
            s ="%3.8f%5s%1.15f\n" %(times[i]," ",Hamils[i])
            file1.write(s)
"""    
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
"""