import csv
from numpy.linalg import norm
from scipy import *
from pylab import plot, show, legend,xlim,ylim,savefig,title,xlabel,ylabel,clf, loglog
import os

"""
wdirbase = "../../../../data/raw/DBsolChrisNoLimiter/o2/"
sdirbase = "../../../../data/postprocessing/DBsolChrisNoLimiter/dx=0.012/"


thetas = [0.0]
#thetas = [1.0]
#for i in range(0,11):
for theta in thetas:
    wdir = wdirbase + "theta="+str(theta) + "/"    
    sdir = sdirbase + "theta="+str(float(theta)) + "/"   
    
    if not os.path.exists(sdir):
        os.makedirs(sdir)
    
             
    s = wdir + "saveoutputtslast.txt"
    with open(s,'r') as file1:
        readfile = csv.reader(file1, delimiter = ',', quotechar='|', quoting=csv.QUOTE_MINIMAL)
            
        h = []
        u = []
        he = []
        ue = []
        x = []
        j = -1
        for row in readfile:       
            if (j >= 0):
                dx =float(row[0])
                dt =float(row[1])
                t =float(row[2])
                Evalf =float(row[4])
                x.append(float(row[3]))
                h.append(float(row[5]))
                u.append(float(row[7]))
                he.append(float(row[8]))
                ue.append(float(row[9]))
                    
            j = j + 1
        x = array(x)
        u = array(u)
        h = array(h)
        ue = array(ue)
        he = array(he)    

   
    n = len(x)
    s = sdir + "NumericalXAndH.dat"
    with open(s,'w') as file1:
        for i in range(n):
            s ="%3.8f%5s%1.15f\n" %(x[i]," ",h[i])
            file1.write(s)
            
    s = sdir + "NumericalXAndU.dat"
    with open(s,'w') as file1:
        for i in range(n):
            s ="%3.8f%5s%1.15f\n" %(x[i]," ",u[i])
            file1.write(s)
            
    s = sdir + "AnalyticalXAndH.dat"
    with open(s,'w') as file1:
        for i in range(n):
            s ="%3.8f%5s%1.15f\n" %(x[i]," ",he[i])
            file1.write(s)
            
    s = sdir + "AnalyticalXAndU.dat"
    with open(s,'w') as file1:
        for i in range(n):
            s ="%3.8f%5s%1.15f\n" %(x[i]," ",ue[i])
            file1.write(s)

s = wdirbase + "savenorms.txt"
with open(s,'r') as file1:
    readfile = csv.reader(file1, delimiter = ',', quotechar='|', quoting=csv.QUOTE_MINIMAL)
        
    dxs = []
    thetas = []
    L1hs = []
    L1us = []
    H1s = []
    j = -1
    for row in readfile:       
        if (j >= 0):
            dxs.append(float(row[0]))
            thetas.append(float(row[1]))
            L1hs.append(float(row[2]))
            L1us.append(float(row[3]))
            H1s.append(abs(float(row[4])))
                
        j = j + 1    
            
s = sdirbase + "ThetaAndHamiltonianConservation.dat"
n = len(thetas)
with open(s,'w') as file1:
    for i in range(n):
        s ="%3.8f%5s%1.15f\n" %(thetas[i]," ",H1s[i])
        file1.write(s)  
        
s = sdirbase + "ThetaAndL1h.dat"
with open(s,'w') as file1:
    for i in range(n):
        s ="%3.8f%5s%1.15f\n" %(thetas[i]," ",L1hs[i])
        file1.write(s) 

s = sdirbase + "ThetaAndL1u.dat"
with open(s,'w') as file1:
    for i in range(n):
        s ="%3.8f%5s%1.15f\n" %(thetas[i]," ",L1us[i])
        file1.write(s) 
"""

wdirbase = "../../../../data/raw/DBNoLimiter/o2/"
sdirbase = "../../../../data/postprocessing/CDBNoLimiter/dx=0.012/"


thetas = [0]
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
