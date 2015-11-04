# -*- coding: utf-8 -*-
"""
Created on Fri May 15 11:23:06 2015

@author: jordan
"""
from scipy import *
from pylab import plot, show, legend,xlim,ylim,savefig,title,xlabel,ylabel,clf, loglog
import csv
from numpy.linalg import norm

def makevar(sx,ex,dx): 
    x = arange(sx, ex, dx)
    return x 

def compareo1(lu,hu,lx,hx,ldx,hdx):
    
    n = len(lx)
    m = len(hx)
    
    ilu = zeros(m)
    
    #boundary case i = 0, should both have same start
    
    ilu[0] = lu[0]
    ilu[-1] = lu[-1]
    
    step = int(ldx / hdx)
    
    k = 1
    for i in range(1,n):
        for j in range(k,k+step):
            ilu[j] = lu[i]
        k = k + step
                          
    return ilu
    

def compareo2(lu,hu,lx,hx,ldx,hdx):
    
    n = len(lx)
    m = len(hx)
    
    ilu = zeros(m)
    ildx = 1.0 / ldx
    
    #boundary case i = 0, should both have same start
    
    ilu[0] = lu[0]
    ilu[-1] = lu[-1]
    
    step = int(ldx / hdx)
    
    k = 1
    for i in range(1,n):
        for j in range(k,k+step):
            diff = hx[j] - lx[i]
            if (i == n-1):
                deriv = 0.5*ildx*(lu[i] - lu[i-1])
            else:
                deriv = 0.5*ildx*(lu[i+1] - lu[i-1])
            ilu[j] = lu[i] + diff*deriv
        k = k + step
                          
    return ilu
                
def dbinit(hl,hf,hc,x):
    n = len(x)
    u = zeros(n)
    
    for i in range(n):
        if (x[i] <= hc):
            u[i] = hl
        else:
            u[i] = hf
    return u

def readdb(wdir,filen):
    s = wdir + filen
    with open(s,'r') as file1:
         readfile = csv.reader(file1, delimiter = ',', quotechar='|', quoting=csv.QUOTE_MINIMAL)
        
         h = []
         u = []
         x = []
         G = []
         j = -1
         for row in readfile:       
             if (j >= 0):
                dx = float(row[0])
                dt = float(row[1])
                t = float(row[2])
                x.append(float(row[3]))
                h.append(float(row[4]))
                G.append(float(row[5])) 
                u.append(float(row[6]))   
                
             j = j + 1
    return dx,dt,t,array(x), array(h), array(u)

wdir = "../../../../data/Joe1/db/o2/"
sdir = "../../../results/bdo2c/"
"""
filen = "0.txt"
ldx,ldt,lt,lx,lh,lu = readdb(wdir,filen)

filen = "5.txt"
hdx,hdt,ht,hx,hh,hu = readdb(wdir,filen)

ilh =  compareo2(lh,hh,lx,hx,ldx,hdx)

plot(hx, hh,'b')
plot(lx,lh, 'r')
plot(hx, ilh, '--k')

#plot(hx, hh - ilh)
normh = norm(hh - ilh, ord = 1) / norm(hh, ord = 1)
"""
files = range(16)
m = len(files)
normhs = zeros(m)
normvs = zeros(m)
ilhs = []
ilus = []
dxs = zeros(m)
s = sdir + "norms2.txt"
with open(s,'a') as file2:
    writefile2 = csv.writer(file2, delimiter = ',', quotechar='|', quoting=csv.QUOTE_MINIMAL)
            
    writefile2.writerow(["dxs" ,"height norm","velocity norm"])   
filen = "15.txt"
hdx,hdt,ht,hx,hh,hu = readdb(wdir,filen)
for i in range(m-7,m):
    filen = str(files[i]) + ".txt"
    ldx,ldt,lt,lx,lh,lu = readdb(wdir,filen)
    dxs[i] = ldx
    ilh =  compareo2(lh,hh,lx,hx,ldx,hdx)
    ilhs.append(ilh)
    ilu =  compareo2(lu,hu,lx,hx,ldx,hdx)
    ilus.append(ilu)
    s = str(ldx)
    plot(lx,lh,label=s)
    print(i,ldx,ldt,ldt/ldx)
    
    normhs[i] = norm(hh - ilh, ord = 1) / norm(hh, ord = 1)
    normvs[i] = norm(hu - ilu, ord = 1) / norm(hu, ord = 1)
    s = sdir + "norms2.txt"
    with open(s,'a') as file2:
        writefile2 = csv.writer(file2, delimiter = ',', quotechar='|', quoting=csv.QUOTE_MINIMAL)
            
        writefile2.writerow([str(ldx) ,str(normhs[i]),str(normvs[i])])  
    legend()
