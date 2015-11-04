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
    
    for i in range(n):
        for j in range(m):
            #print (i,j)
            diff = lx[i] - hx[j] 
            if((diff > - 0.5*ldx) and (diff < 0.5*ldx)): 
                #then we are inside the lower cell with middle lx[i]
                ilu[j] = lu[i]
            elif(diff == - 0.5*ldx):
                if (i == 0):
                    ilu[j] = lu[i]
                else:
                    if(abs(lu[i] - hu[j]) >abs(lu[i-1] - hu[j])):
                        ilu[j] = lu[i-1]
                    else:
                        ilu[j] = lu[i]
            elif(diff == 0.5*ldx):
                if (i == n-1):
                    ilu[j] = lu[i]
                else:
                    if(abs(lu[i] - hu[j]) >abs(lu[i+1] - hu[j])):
                        ilu[j] = lu[i+1]
                    else:
                        ilu[j] = lu[i]
                    
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

wdir = "../data/Cserre/dambreak/o1r/"
sdir = "../data/Cserre/dambreak/o1r/"

#filen = "28.txt"
#ldx,ldt,lt,lx,lh,lu = readdb(wdir,filen)

#filen = "27.txt"
#hdx,hdt,ht,hx,hh,hu = readdb(wdir,filen)

#ilh =  compareo1(lh,hh,lx,hx,ldx,hdx)

#plot(hx, hh,'b')
#plot(lx,lh, 'r')
#plot(hx, ilh, '--k')

#plot(hx, hh - ilh)
#normh = norm(hh - ilh, ord = 1) / norm(hh, ord = 1)

files = range(22)
m = len(files)
normhs = zeros(m)
normvs = zeros(m)
ilhs = []
ilus = []
dxs = zeros(m)
s = wdir + "norms.txt"
with open(s,'a') as file2:
    writefile2 = csv.writer(file2, delimiter = ',', quotechar='|', quoting=csv.QUOTE_MINIMAL)
            
    writefile2.writerow(["dxs" ,"height norm","velocity norm"])   
filen = "36.txt"
hdx,hdt,ht,hx,hh,hu = readdb("../data/Cserre/dambreak/o1/",filen)
for i in range(m):
    filen = str(files[i]) + ".txt"
    ldx,ldt,lt,lx,lh,lu = readdb(wdir,filen)
    dxs[i] = ldx
    ilh =  compareo1(lh,hh,lx,hx,ldx,hdx)
    ilhs.append(ilh)
    ilu =  compareo1(lu,hu,lx,hx,ldx,hdx)
    ilus.append(ilu)
    normhs[i] = norm(hh - ilh, ord = 1) / norm(hh, ord = 1)
    normvs[i] = norm(hu - ilu, ord = 1) / norm(hu, ord = 1)
    s = wdir + "norms.txt"
    with open(s,'a') as file2:
        writefile2 = csv.writer(file2, delimiter = ',', quotechar='|', quoting=csv.QUOTE_MINIMAL)
            
        writefile2.writerow([str(ldx) ,str(normhs[i]),str(normvs[i])])         
