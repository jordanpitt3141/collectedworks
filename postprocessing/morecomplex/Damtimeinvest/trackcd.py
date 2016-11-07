# -*- coding: utf-8 -*-
"""
Created on Wed Sep 23 16:53:37 2015

@author: jordan
"""

import csv
from numpy.linalg import norm
from scipy import *
from pylab import plot, show, legend,xlim,ylim,savefig,title,xlabel,ylabel,clf, loglog, xticks,yticks
from matplotlib2tikz import save as tikz_save
import os
from subprocess import call


diffs = [10.0]
sdirbase = "../../../../data/postprocessing/uhcomp/DB/o3/10/"
wdirbase = "../../../../data/raw/trackleadsola10new/o3/"

def centreddiff(x,q,dx):
    idx = 1.0 / dx
    n = len(q)
    dq = zeros(n)
    for i in range(1, n-1):
        dq[i] =0.5*idx*(q[i+1] - q[i-1])
        
    dq[0] =0.5*idx*(q[1] - q[0])
    
    dq[n-1] =0.5*idx*(q[n-1] - q[n-2])
    
    return dq
    
def findzeros(q):
    n = len(q)
    zeros = []
    for i in range(1,n-1):
        if(q[i+1]*q[i-1] <= 0):
            zeropt = i 
            zeros.append(zeropt)
    return zeros
        


hps = [] 
ups = []
ts = []   
#for i in range(5120,1024000,5120):
for i in range(0,1):
    
    i = 1024000
    

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
                x.append(float(row[4]))
                h.append(float(row[5]))
                u.append(float(row[7]))
             j = j + 1
                     
    u2 = 1.074975
    h2 = 1.36898
    
    dh = centreddiff(x,h,dx)
    du = centreddiff(x,u,dx)
    
    dhz = findzeros(dh)
    duz = findzeros(du)
    
    x2 = 500 + u2*t
    x2i = int(x2/dx)
    hp = h[x2i]
    up = u[x2i]
    hps.append(hp)
    ups.append(up)
    ts.append(t)
    
h2const = h2*ones(len(ts))
u2const = u2*ones(len(ts))