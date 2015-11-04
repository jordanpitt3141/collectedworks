# -*- coding: utf-8 -*-
"""
Created on Wed Sep 23 16:53:37 2015

@author: jordan
"""

import csv
from numpy.linalg import norm
from scipy import *
from pylab import plot, show, legend,xlim,ylim,savefig,title,xlabel,ylabel,clf, loglog
from matplotlib2tikz import save as tikz_save

diffw = "11"
wdirord = "FDcent"

#nums = [0,3,6,9,12]
#nums = [2,5,8,11,14]
#nums = [1,4,7,10,13]
nums = range(9,13)
for k in nums:
    wdir = "../../../../../data/Joesmooth/bigsmooth/"  +wdirord +"/" + str(k) + "/" + diffw + "/" 
    sdir = "../../../../results/smoothdb/1/"
         
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
                u.append(float(row[5])) #5 for FDcent, 6 otherwise 
                beta = float(row[6]) #could be 8 as well for o2
                
             j = j + 1
         x = array(x)
         h = array(h)     
  
    s = str(dx) 
    plot(x,h ,label=s)
    ylim([-0.1,2])
    xlim([0,1000])
    s = "Dam Break: " + wdirord + " diff = " + str(beta)
    title(s)
    xlabel("x (m)")
    ylabel("h (m)")
    legend()
        
    """
    plot(x,u ,'-b')
    ylim([-0.1,2])
    xlim([0,1000])
    s = "Dam Break: " + wdirord + " dx = " + str(dx)
    title(s)
    xlabel("x (m)")
    ylabel("u (m/s)")
    """
    
s = sdir + wdirord +".tikz" 
tikz_save(s);      
clf();