import csv
from numpy.linalg import norm
from scipy import *
from pylab import plot, show, legend,xlim,ylim,savefig,title,xlabel,ylabel,clf, loglog
import os

wdir = "../../../../data/raw/SolColAGN/o3/"
sdir = "../../../../data/postprocessing/SolColAGN/o3/"


if not os.path.exists(sdir):
    os.makedirs(sdir)

s = wdir + "saveoutputtslast.txt"
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
            Evaln = float(row[3])
            x.append(float(row[4]))
            h.append(float(row[5]))
            u.append(float(row[7]))
               
        j = j + 1

s = wdir + "saveoutputts1.txt"
with open(s,'r') as file1:
    readfile = csv.reader(file1, delimiter = ',', quotechar='|', quoting=csv.QUOTE_MINIMAL)

    j = -1
    for row in readfile:       
        if (j >= 0):
            Evali = float(row[3])
               
        j = j + 1

n = len(x)
s = sdir + "h.dat"
with open(s,'w') as file1:
    for i in range(n):
        s ="%3.8f%5s%1.15f\n" %(x[i]," ",h[i])
        file1.write(s)
s = sdir + "u.dat"
with open(s,'w') as file2:
    for i in range(n):
        s ="%3.8f%5s%1.15f\n" %(x[i]," ",u[i])
        file2.write(s)

s = sdir + "Energy.dat"
with open(s,'w') as file1:
    s ="%11s%5s%16s\n" %("dx"," ","H(0)")
    file1.write(s)  
    s ="%3.8f%5s%1.15f\n" %(dx," ",Evali)
    file1.write(s)      
    s ="%11s%5s%16s\n" %("dx"," ","H(50s)")
    file1.write(s)  
    s ="%3.8f%5s%1.15f\n" %(dx," ",Evaln)
    file1.write(s)    
    s ="%11s%5s%16s\n" %("dx"," ","(H(0) - H(50s))/H(0)")
    file1.write(s)  
    s ="%3.8f%5s%1.15f\n" %(dx," ",(Evali - Evaln)/ Evali)
    file1.write(s)      


