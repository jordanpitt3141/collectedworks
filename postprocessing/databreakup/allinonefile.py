# -*- coding: utf-8 -*-
"""
Created on Tue Jul  5 10:47:33 2016

@author: jordan
"""

import csv
from numpy.linalg import norm
from scipy import *
from pylab import plot, show, legend,xlim,ylim,savefig,title,xlabel,ylabel,clf, loglog
from numpy import ones
import os

wdir = "../../../data/raw/NEWdata/DSWalphalongt/o3/10/0/"

sdirb = "../../../data/raw/DSWalphalongfix/o3/10/"

fils = os.listdir(wdir)

for fil in fils:
    print(fil)
    s = wdir + fil
    with open(s,'r') as file1:
        readfile = csv.reader(file1, delimiter = ',', quotechar='|', quoting=csv.QUOTE_MINIMAL)
            
        hs = []
        Gs = []
        us = []
        xs = []
        dxs = []
        dts = []
        ts = []
        Es = []
        diffs = []
        dx = 0.0
        dt = 0.0
        t = 0.0
        E = 0
        x = []
        h = []
        G = []
        u = []
        diff = 0
        for row in readfile:
            if row[0] == 'dx':
                
                dxs.append(dx)
                dts.append(dt)
                ts.append(t)
                Es.append(E)
                xs.append(x)
                hs.append(h)
                Gs.append(G)
                us.append(u)
                diffs.append(diff)
                
                x = []
                h = []
                G = []
                u = []
                
            else:            
                dx =float(row[0])
                dt =float(row[1])
                t =float(row[2])
                E = float(row[3])
                x.append(float(row[4]))
                h.append(float(row[5]))
                G.append(float(row[6]))
                u.append(float(row[7]))
                diff = float(row[8])
    
    dxs.append(dx)
    dts.append(dt)
    ts.append(t)
    Es.append(E)
    xs.append(x)
    hs.append(h)
    Gs.append(G)
    us.append(u)
    diffs.append(diff)  
    
    for k in range(1,len(dxs)):
        sdir = sdirb + str(diffs[k]) + "/"
        
        if not os.path.exists(sdir):
            os.makedirs(sdir) 
        s = sdir + fil
        with open(s,'a') as file2:
             writefile2 = csv.writer(file2, delimiter = ',', quotechar='|', quoting=csv.QUOTE_MINIMAL)
            
             writefile2.writerow(['dx' ,'dt','time','Eval','cell midpoint', 'height(m)', 'G' , 'u(m/s)',"diffuse"])        
                           
             for j in range(len(xs[1])):
                 writefile2.writerow([str(dxs[k]),str(dts[k]),str(ts[k]),str(Es[k]), str(xs[k][j]), str(hs[k][j]) , str(Gs[k][j]) , str(us[k][j]), str(diffs[k])])   