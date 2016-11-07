# -*- coding: utf-8 -*-
"""
Created on Thu Sep 11 16:26:29 2014

@author: Jordan
"""

import csv
from numpy.linalg import norm
from scipy import *
from pylab import plot, show, legend,xlim,ylim,savefig,title,xlabel,ylabel,clf, loglog

betanum = "10"
#wdirord = "o3"

#wdir = "../../../data/Joe/alldb/" +wdirord+"/" 
wdir = "../../../data/t/3/" 
sdir = "../../results/show/steve1/correct/"
     
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
            u.append(float(row[5])) 
            #beta = float(row[7])
            
         j = j + 1
     x = array(x)
     h = array(h)     
"""
plot(x,h ,'-b')
ylim([-0.1,2])
xlim([0,1000])
s = "Dam Break: " + wdirord + " dx = " + str(dx)
title(s)
xlabel("x (m)")
ylabel("h (m)")
"""
plot(x,u ,'-b')
ylim([-0.1,2])
xlim([0,1000])
s = "Dam Break: " + wdirord + " dx = " + str(dx)
title(s)
xlabel("x (m)")
ylabel("u (m/s)")