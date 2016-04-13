# -*- coding: utf-8 -*-
"""
Created on Thu Sep 11 16:26:29 2014

@author: Jordan
"""

import csv
from numpy.linalg import norm
from scipy import *
from pylab import plot, show, legend,xlim,ylim,savefig,title,xlabel,ylabel,clf, loglog
from Hamil import *
from os import listdir

def copyarraytoC(a):
    n = len(a)
    b = mallocPy(n)
    for i in range(n):
        writetomem(b,i,a[i])
    return b
    
def copyarrayfromC(a,n):
    b = [0]*n
    for i in range(n):
        b[i] = readfrommem(a,i)
        
    return b

dxw = "3"
wdirord = "o1"
wdatadir = "../../../../../data/raw/bigsmoothtargetted/"  +wdirord +"/"


dxws = listdir(wdatadir)
Evals = []


for dxw in dxws:
    wdatadirn = wdatadir + dxw + "/" 
    nums = listdir(wdatadirn)
    for k in nums:
        wdir = wdatadirn  + str(k) + "/" 
        sdir = "../../../results/dbEnergy/" +wdirord +"/" + dxw + "/" + str(k) + "/"
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
                    dx = float(row[0])
                    dt = float(row[1])
                    t = float(row[2])
                    x.append(float(row[3]))
                    h.append(float(row[4]))
                    u.append(float(row[6])) 
                    beta = float(row[7]) #could be 8 as well for o2
                    
                 j = j + 1
             x = array(x)
             h = array(h)     
      
        startx = x[0]
        endx = x[-1]
        g = 9.81
        n = len(x)
        niBC = 3    
        xbeg = arange(startx - niBC*dx,startx,dx)
        xend = arange(endx + dx,endx + (niBC+1)*dx) 
        hbeg = h[0]*ones(niBC)
        hend = h[-1]*ones(niBC)
        ubeg = u[0]*ones(niBC)
        uend = u[-1]*ones(niBC)
        xbc =  concatenate([xbeg,array(x),xend])
        hbc =  concatenate([hbeg,array(h),hend])
        ubc =  concatenate([ubeg,array(u),uend])
        
        xbc_c = copyarraytoC(xbc)
        hbc_c = copyarraytoC(hbc)
        ubc_c = copyarraytoC(ubc)
        
        Eval = HankEnergyall(xbc_c,hbc_c,ubc_c,g,n + 2*niBC,niBC,dx)
        Evals.append((Eval,wdirord,dxw,beta))