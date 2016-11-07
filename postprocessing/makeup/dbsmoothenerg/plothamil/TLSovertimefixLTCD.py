import csv
from numpy.linalg import norm
from scipy import *
from pylab import plot, show, legend,xlim,ylim,savefig,title,xlabel,ylabel,clf, loglog
from os import listdir
import os

wdirb = "../../../../../data/postprocessing/RLHamil/"
sdir = "../../../../../data/postprocessing/RLHamilOT/"

if not os.path.exists(sdir):
    os.makedirs(sdir) 

times = arange(0,300.5,0.5)
t = []
Hams = []
FTs = []
STs = []
TTs = []
Hamis = []
FTis = []
STis = []
TTis = []

for ts in times:
    wdir = wdirb + str(ts) + "secs/"         
    s = wdir + "relperc.txt"
    with open(s,'r') as file1:
        readfile = csv.reader(file1, delimiter = ',', quotechar='|', quoting=csv.QUOTE_MINIMAL)
        j = -1
        for row in readfile:
            if (j >= 0):
                dx = float(row[0])
                dt =float(row[1])
                t.append(float(row[2]))
                beta = float(row[3])
                Hams.append(float(row[4]))
                FTs.append(float(row[5]))
                STs.append(float(row[6]))
                TTs.append(float(row[7]))
                Hamis.append(float(row[8]))
                FTis.append(float(row[9]))
                STis.append(float(row[10]))
                TTis.append(float(row[11]))
            j = j + 1

n = len(t)
s = sdir + "FT.dat"
with open(s,'w') as file1:
    s ="%3.8f%5s%1.15f\n" %(0.0," ",FTis[0])
    file1.write(s)
    for i in range(n):
        s ="%3.8f%5s%1.15f\n" %(t[i]," ",FTs[i])
        file1.write(s)
s = sdir + "ST.dat"
with open(s,'w') as file1:
    s ="%3.8f%5s%1.15f\n" %(0.0," ",STis[0])
    file1.write(s)
    for i in range(n):
        s ="%3.8f%5s%1.15f\n" %(t[i]," ",STs[i])
        file1.write(s)
s = sdir + "TT.dat"
with open(s,'w') as file1:
    s ="%3.8f%5s%1.15f\n" %(0.0," ",TTis[0])
    file1.write(s)
    for i in range(n):
        s ="%3.8f%5s%1.15f\n" %(t[i]," ",TTs[i])
        file1.write(s)
       
         


  