# -*- coding: utf-8 -*-
"""
Created on Wed Sep 23 16:53:37 2015

@author: jordan
"""

import csv
from numpy.linalg import norm
from scipy import *

from pylab import plot, show, legend,xlim,ylim,savefig,title,xlabel,ylabel,clf, loglog
import os
from subprocess import call


diffs = [10.0]
sdirbase = "../../../../data/postprocessing/RLTimeAplus/9/"
wdirbase = "../../../../data/raw/longcontactdiscdx9diff10fileio/o3/9/0/"
#wdirbase = "../../../../data/raw/DBASPECTRAT/o3/10/10/8.0/"

def centreddiff(x,q,dx):
    idx = 1.0 / dx
    n = len(q)
    dq = zeros(n)
    for i in range(1, n-1):
        dq[i] =0.5*idx*(q[i+1] - q[i-1])
        
    dq[0] =0.5*idx*(q[1] - q[0])
    
    dq[n-1] =0.5*idx*(q[n-1] - q[n-2])
    
    return dq
    
def findzeros(q,Q,x):
    n = len(q)
    qxs = []
    qvs = []
    Qvs = []
    for i in range(1,n):
        if(q[i]*q[i-1] <= 0 and q[i] < q[i-1] ):
            qx = 0.5*(x[i] + x[i-1])
            qv = 0.5*(q[i] + q[i-1])
            Qv = 0.5*(Q[i] + Q[i-1])
            qxs.append(qx)
            qvs.append(qv)
            Qvs.append(Qv)
    return qxs,qvs,Qvs
    
def findleadsol(dqx,qv,x,h,dx):
    n = len(qv)
    nqv = []
    ndqx = []
    for i in range(1,n):
            if(qv[i] > 1.1 and qv[i] < 1.79):
                nqv.append(qv[i])
                ndqx.append(dqx[i])
    if (len(nqv) == 0):
        nqv = [0]
        ndqx = [500]
    leadsolh = nqv[-1]
    leadsolx = ndqx[-1]
    tol = 1
    xlow = int((leadsolx - tol + 900)/dx)
    xhi = int((leadsolx + tol + 900)/dx)
   
    leadsolh = 0
    for i in range(xlow,xhi):
        
        if(h[i] > leadsolh):
            leadsolh = h[i] 
            leadsolx = x[i]
        
    return leadsolh, leadsolx


ts = [] 
aplus = []


if not os.path.exists(sdirbase):
    os.makedirs(sdirbase)
for i in range(0,300*10,5):
#for i in range(0,1):
        
    s = wdirbase + "out" + str(i) + ".txt"
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
             j = j + 1
                     
    u2 = 1.074975
    h2 = 1.36898
    
    dh = centreddiff(x,h,dx)    
    dhx,dhv,hv = findzeros(dh,h,x)
    
    solqv , solqx = findleadsol(dhx,hv,x,h,dx)
    
    aplus.append(solqv)
    ts.append(t)
    
    #plot(x,h)
    #xlim([dhx[-1] - 10 , dhx[-1] + 10])
    
n = len(ts)
s = sdirbase + "aplus.dat"
with open(s,'w') as file1:
    for i in range(n):
        s ="%3.8f%5s%1.15f\n" %(ts[i]," ",aplus[i])
        file1.write(s)


    
    
