# -*- coding: utf-8 -*-
"""
Created on Wed Sep 23 16:53:37 2015

@author: jordan
"""

import csv
from numpy.linalg import norm
from scipy import *
from pylab import plot, show, legend,xlim,ylim,savefig,title,xlabel,ylabel,clf, loglog
#from matplotlib2tikz import save as tikz_save
import os

dxw = "5"
wdirord = "o1"

wdirords = ["o3","FDcent","grim","o2","o1"]
ordsubtup = [[6,7],[5,6],[5,6], [6,8], [6,7]]
removeint = [[515,545],[515,545],[515,545],[515,545],[515,545]]

dxws = range(3,11)

numdiff = range(21)

#wdirords = ["o3"]
#ordsubtup = [[6,7]]
#removeint = [[515,545]]

#dxws = [8]
#numdiff = [13,14,15]

for ip in range(len(wdirords)):
    for jp in dxws:
        
        hs = []
        us = []
        xs = []
        ts = []
        dxs = []
        dts = []
        diffs = []
        normhs =[]
        normus =[]
        
        dxw = str(jp)
        wdirord = wdirords[ip]
        
        nums = numdiff
                    
        for k in nums:
            wdir = "../../../../../data/Joesmooth/bigsmooth/"  +wdirord +"/" + dxw + "/" + str(k) + "/"
            sdir = "../../../../results/smoothdb/1/1dxmdiffcomreal/" + wdirord + "/" +dxw+ "/"
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
                        u.append(float(row[ordsubtup[ip][0]])) #5 for FDcent, 6 otherwise 
                        beta = float(row[ordsubtup[ip][1]]) #could be 8 as well for o2
                    j = j + 1
                    
                zbeg = int(removeint[ip][0]/dx)
                zend = int(removeint[ip][1]/dx)
                h = h[0:zbeg] + zeros(zend-zbeg).tolist()  + h[zend:]
                u = u[0:zbeg] + zeros(zend-zbeg).tolist()  + u[zend:]
                x = array(x)
                h = array(h)
                u = array(u)
                hs.append(h)
                us.append(u)
                xs.append(x)
                dxs.append(dx)
                dts.append(dt)
                ts.append(t)
                diffs.append(beta)
        

        s = sdir + "norms.txt"
        with open(s,'a') as file2:
            writefile2 = csv.writer(file2, delimiter = ',', quotechar='|', quoting=csv.QUOTE_MINIMAL)
                    
            writefile2.writerow(["dx","time","diff" ,"height norm","velocity norm"])           
        # Comparison
        for i in numdiff:
            j = i - numdiff[0]
            #plot(xs[j],hs[j])
            normh = norm(hs[j] - hs[-1],ord=1) / norm(hs[-1],ord=1)
            normu = norm(us[j] - us[-1],ord=1) / norm(us[-1],ord=1)
            
            s = sdir + "norms.txt"
            with open(s,'a') as file1:
                writefile1 = csv.writer(file1, delimiter = ',', quotechar='|', quoting=csv.QUOTE_MINIMAL)
                    
                writefile1.writerow([str(dxs[j]),str(ts[j]),str(diffs[j]) ,str(normh),str(normu)]) 
            
            normus.append(normu)
            normhs.append(normh)
                
                
